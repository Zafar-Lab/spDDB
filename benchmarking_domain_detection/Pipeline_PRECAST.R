library(PRECAST)
library(Seurat)
library(scater)
library(anndata)
set.seed(101)

ST <- read_h5ad("ST.h5ad")
ST$var_names_make_unique()

count <- t(ST$X)
meta_data <- data.frame(row=ST$obs$array_row, col=ST$obs$array_col)
row.names(meta_data) <- colnames(count)

seu <- CreateSeuratObject(counts=count, meta.data = meta_data)

meta_data <- seu@meta.data
all(c("row", "col") %in% colnames(meta_data))  ## the names are correct!

head(meta_data[, c("row", "col")])

preobj <- CreatePRECASTObject(seuList = list(seu), selectGenesMethod = "HVGs", gene.number = 2000)

PRECASTObj <- AddAdjList(preobj, platform = "Visium")

PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, coreNum = 1, maxIter = 30, verbose = TRUE)

q <- length(unique(ST$obs$ground_truth))
PRECASTObj <- PRECAST(PRECASTObj, K = q)

#Put the reults into a Seurat object seuInt.
seuInt <- PRECASTObj@seulist[[1]]
seuInt@meta.data$cluster <- factor(unlist(PRECASTObj@resList[[1]]$cluster[[1]]))
seuInt@meta.data$batch <- 1
seuInt <- Add_embed(PRECASTObj@resList[[1]]$hZ[[1]], seuInt, embed_name = "PRECAST")
posList <- lapply(PRECASTObj@seulist, function(x) cbind(x$row, x$col))
seuInt <- Add_embed(posList[[1]], seuInt, embed_name = "position")
Idents(seuInt) <- factor(seuInt@meta.data$cluster)

#saving results
final <- seuInt@meta.data
key <- seuInt$row
layer_guess <- ST$obs$ground_truth[key]
final['annotation'] <- layer_guess

write.csv(final, "/home/keerthana/Documents/Spatial/spatial_domain/final_outputs/mouse_breast_cancer/output_PRECAST.csv")
