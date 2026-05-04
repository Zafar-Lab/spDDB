library(SpatialDecon)
library(anndata)
library(Seurat)
library(scater)

h5adpath <- "./simulated_gene_expression_average.h5ad" 
dlpfc <- read_h5ad(h5adpath)
st_count <- t(dlpfc$X)

st_count <- round(st_count)

andersson_g1 = CreateSeuratObject(counts = st_count, assay="Spatial")

pos <- dlpfc$obsm$spatial
andersson_g1@meta.data$x = pos[,1]
andersson_g1@meta.data$y = pos[,2]

SC <- read_h5ad("./SC.h5ad")
mtx <- t(SC$X)

cell_type <- SC$obs$cell_type
cell_type <- factor(cell_type)
CellID <- factor(mtx@Dimnames[[2]])
sc_annmeta<- data.frame(CellID = CellID, LabeledCellType= cell_type)

custom_mtx <- create_profile_matrix(mtx =mtx,            # cell x gene count matrix
                                    cellAnnots = sc_annmeta,  # cell annotations with cell type and cell name as columns 
                                    cellTypeCol = "LabeledCellType",  # column containing cell type
                                    cellNameCol = "CellID",           # column containing cell ID/name
                                    matrixName = "custom_mini_colon", # name of final profile matrix
                                    outDir = NULL,                    # path to desired output directory, set to NULL if matrix should not be written
                                    normalize = FALSE,                # Should data be normalized? 
                                    minCellNum = 5,                   # minimum number of cells of one type needed to create profile, exclusive
                                    minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
                                    scalingFactor = 5,                # what should all values be multiplied by for final matrix
                                    discardCellTypes = TRUE)  

res = runspatialdecon(object = andersson_g1,
                      # bg = 0.01,  
                      X = custom_mtx,
                      align_genes = TRUE)

path <- paste("output.csv")

final <- t(res$prop_of_all)
final <- final[ , sort(colnames(final))]
write.csv(final , path)


