library("DR.SC")
library(Seurat)
library(anndata)
library("scater")

ST <- read_h5ad("./ST.h5ad")
ST$var_names_make_unique()

count <- t(ST$X)
meta_data <- data.frame(row=ST$obsm$spatial[,1], col=ST$obsm$spatial[,2])
row.names(meta_data) <- colnames(count)

seu <- CreateSeuratObject(counts=count, meta.data = meta_data)


#######
seuf <- NormalizeData(seu, verbose = F)
# choose 500 highly variable features
seuf <- FindVariableFeatures(seuf, nfeatures = 500, verbose = F)

q <- length(unique(ST$obs$annotation))
seuf <- DR.SC(seuf, K= q, platform = 'Visium', verbose=F)
result <- data.frame('cluster' = seuf$spatial.drsc.cluster , row.names = colnames(seuf))

write.csv(result ,"out.csv")