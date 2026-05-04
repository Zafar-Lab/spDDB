library(Seurat)
library(scater)
library(iSC.MEB)
library(SingleCellExperiment)
library(anndata)

path <- "./ST.h5ad"
ST <- read_h5ad(path)

#Number of clusters
q <- length(unique(ST$obs$ground_truth))
ST$var_names_make_unique()
counts <- t(ST$X)

x<- ST$obs$array_row
y<- ST$obs$array_col
spatial_location_anndata<- data.frame(row =x, col=y)
row.names(spatial_location_anndata)<-colnames(counts)
spatial_location_anndata <- as.matrix(spatial_location_anndata)

ST <- SingleCellExperiment(assays = list(counts = t(ST$X)) , colData = spatial_location_anndata , rowData = ST$var)
ST <- logNormCounts(ST)

ST.seurat <- as.Seurat(ST )
ST.seurat <- RenameAssays(ST.seurat , assay.name = "originalexp" , new.assay.name =  "RNA")

seulist <- list(i = ST.seurat)

iSCMEBObj <- CreateiSCMEBObject(seuList = seulist, verbose = TRUE, premin.spots = 0, postmin.spots = 0)

iSCMEBObj <- CreateNeighbors(iSCMEBObj, platform = "Visium")

iSCMEBObj <- runPCA(iSCMEBObj, npcs = 15, pca.method = "APCA")

iSCMEBObj <- SetModelParameters(iSCMEBObj, verbose = TRUE)

iSCMEBObj <- iSCMEB(iSCMEBObj, K = q)

layer_pred <- idents(iSCMEBObj)
key <- iSCMEBObj@seulist$i$row
result <- data.frame( cluster = layer_pred[[1]])
row.names(result) <- names(key)
name <- "./output_ISC_MEB.csv"
write.csv(result , file = name)

