library(Giotto)
library(anndata)
library(Matrix)
library(Seurat)
library(scran)
library(reticulate)
args<-commandArgs(trailingOnly = TRUE)
sc_path <- "/data/narein/ST_deconvolution_Tanush/Data/MERFISH_Ileum/SC/SC.h5ad"
st_path <- "/data/narein/ST_deconvolution_Tanush/Data/MERFISH_Ileum/ST/ST.h5ad"
output_path <- "/data/narein/ST_deconvolution_Tanush/Data/MERFISH_Ileum/final_outputs/output_SpatialDWLS.csv"
celltype_key <- "cell_type"
sc_path <-args[1]
st_path <- args[2]
celltype_key <- args[3]
output_path <- args[4]
python_path <- "/home/narein/.virtualenvs/r-reticulate/bin/python"

st<-read_h5ad(st_path)
new_x <- st$obsm$spatial[,1]
new_y <- st$obsm$spatial[,2]
st$var_names <- make.unique(st$var_names)
#if ("array_row" %in% colnames(st$obs)) {
#  new_x <- st$obs['array_row']
#  new_y <- st$obs['array_col']
#} else if("new_x" %in% colnames(st$obs)){
#  new_x <- st$obs['new_x']
#  new_y <- st$obs['new_y']
#} else{
#  new_x <- st$obsm$spatial[,1]
#  new_y <- st$obsm$spatial[,2]
#}
coords_df <- data.frame(cbind(new_x,new_y))
st<- as.matrix(st)
st<- t(st)
st <- round(st)
x <- colSums(st) > 0
st <- st[,x]
coords_df <- coords_df[x,]
instrs = createGiottoInstructions(python_path = python_path)
st_giotto <- createGiottoObject(raw_exprs = st,spatial_locs = coords_df,instructions = instrs)
rm(st)
st_giotto <- normalizeGiotto(gobject = st_giotto)
#st_giotto <- calculateHVF(gobject = st_giotto,expression_values = 'normalized',return_plot = FALSE)
gene_metadata = fDataDT(st_giotto)
hvggenes = gene_metadata$feat_ID
#hvggenes = gene_metadata[hvf == 'yes']$gene_ID
st_giotto <- runPCA(gobject = st_giotto, genes_to_use = hvggenes, scale_unit = F)
signPCA(st_giotto,scale_unit = F)
st_giotto <- runUMAP(st_giotto, dimensions_to_use = 1:30)
st_giotto <- createNearestNetwork(gobject = st_giotto, dimensions_to_use = 1:30, k = 15)
st_giotto <- doLeidenCluster(gobject = st_giotto, resolution = 0.3, n_iterations = 1000)

sc<-read_h5ad(sc_path)
sc$var_names <- make.unique(sc$var_names)
celltypes <- sc$obs[celltype_key][,1]
sc <- as.matrix(sc)
sc<- t(sc)
sc_giotto <- createGiottoObject(raw_exprs = sc,instructions = instrs)
rm(sc)
sc_giotto <- normalizeGiotto(gobject = sc_giotto)
sc_giotto@cell_metadata$cell$rna@metaDT$leiden_clus <- celltypes
markers = findMarkers_one_vs_all(gobject = sc_giotto,cluster_column = 'leiden_clus',
                                 method = 'scran',
                                 expression_values = 'normalized')
Sig_scran <- unique(markers$feats[which(markers$ranking <= 20)])
norm_exp<-2^(sc_giotto@expression$cell$rna$normalized@exprMat)-1
id<-sc_giotto@cell_metadata$cell$rna@metaDT$leiden_clus
rm(sc_giotto)
ExprSubset<-norm_exp[Sig_scran,]
Sig_exp<-NULL
for (i in unique(id)){
  Sig_exp<-cbind(Sig_exp,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
}
colnames(Sig_exp)<-unique(id)

result = runDWLSDeconv(gobject = st_giotto, sign_matrix = Sig_exp, return_gobject = FALSE)
rm(st_giotto)
df <- data.frame(result@enrichDT)
rownames(df) <- df$cell_ID
df <- df[-1]
df <-df[,order(colnames(df))]
write.csv(df,output_path)
