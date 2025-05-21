library(anndata)
library(Matrix)
library(Seurat)
library(dplyr)
set.seed(1)

##Read paths
args<-commandArgs(trailingOnly = TRUE)
st_path = args[2]
sc_path <-args[1]
celltype_key <- args[3]
output_path <- args[4]

sc<-read_h5ad(sc_path)
sc <- sc[rowSums(sc$X)!=0]
sc$var_names <- make.unique(sc$var_names)
celltypes <- sc$obs[celltype_key]
sc<-as.matrix(sc)
sc <- t(sc)
sc <- CreateSeuratObject(counts = sc)

st<-read_h5ad(st_path)
st$var_names <- make.unique(st$var_names)
st<-as.matrix(st)
st<-t(st)
#st <- round(st)
st <- st[,colSums(st) > 0]
st<-CreateSeuratObject(st)
st<-SCTransform(st, assay="RNA", verbose=TRUE)
st<-RunPCA(st,assay="SCT",verbose=FALSE)
st <- FindNeighbors(st, reduction = "pca", dims = 1:30)
st <- FindClusters(st,verbose=FALSE)
st<-RunUMAP(st,reduction="pca", dims=1:30)

## Run Deconvolution
sc<-SCTransform(sc,ncells=3000,verbose=FALSE)
sc<-RunPCA(sc,verbose=FALSE)
sc<-RunUMAP(sc,dims=1:30)
anchors <- FindTransferAnchors(reference = sc, query = st, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = celltypes[,1], prediction.assay = TRUE,weight.reduction = st[["pca"]], dims = 1:30)

## Convert output to correct format
res <- t(predictions.assay@data)
column_index_to_delete <- which(colnames(res) == 'max')
res <- res[,-column_index_to_delete]
res <- res[,sort(colnames(res))]
write.csv(res,output_path)
