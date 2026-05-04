library(CARD)
library(scater)
library(anndata)

st_path <- "./simulated_gene_expression_average.h5ad"

ST <- read_h5ad(st_path)
SC <- read_h5ad("SC.h5ad")

sc_adata<- t(SC$X)
st_adata <- t(ST$X)
st_adata <- round(st_adata)

cell_type <- SC$obs$cell_type
cell_type <- factor(cell_type)
cellID <- factor(colnames(sc_adata))
sc_annmeta<- data.frame(cellID = cellID, cell_type= cell_type)
sc_annmeta$sampleInfo = "sample1"
rownames(sc_annmeta) <- colnames(sc_adata)

pos <- ST$obsm$spatial

row.names(pos)<-colnames(st_adata)
colnames(pos) <- c('x','y')
spatial_location_anndata <- as.data.frame(pos)
row.names(spatial_location_anndata)<-colnames(st_adata)

CARD_obj = createCARDObject(
  sc_count = sc_adata,
  sc_meta = sc_annmeta,
  spatial_count = st_adata,
  spatial_location = spatial_location_anndata,
  ct.varname = "cell_type",
  ct.select = unique(sc_annmeta$cell_type),
  sample.varname = "sampleInfo",
  minCountGene = 40,
  minCountSpot = 2)

CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

#saving results
path <- "output.csv"
final <- CARD_obj@Proportion_CARD
final <- final[ , sort(colnames(final))]
write.csv(final,file=path)