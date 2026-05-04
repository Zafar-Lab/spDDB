# remotes::install_github("prabhakarlab/Banksy")
library(Banksy)
library(scater)
library(anndata)
library(SpatialExperiment)
library(cowplot)
set.seed(101)

path <- "./ST.h5ad"
ST <- read_h5ad(path)
st_adata <- t(ST$X)

#Spatial coords
x<- ST$obs$array_row
y<- ST$obs$array_col
spatial_location_anndata<- data.frame(sdimx =x, sdimy=y)
row.names(spatial_location_anndata)<-colnames(st_adata)

spatial_location_anndata <- as.matrix(spatial_location_anndata)

se <- SpatialExperiment(assay = list(counts = st_adata), spatialCoords = spatial_location_anndata)

se <- computeLibraryFactors(se)
aname <- "normcounts"
assay(se, aname) <- normalizeCounts(se, log = FALSE)

#0,8 for MERFISH/osmfish
lambda <- c(0, 0.2)
k_geom <- 18

se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)

se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)

#Resolution varied as per number of clusters
res = 0.3
se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = res)
se <- Banksy::connectClusters(se)

cnames <- colnames(colData(se))
cnames <- cnames[grep("^clust", cnames)]
colData(se) <- cbind(colData(se), spatialCoords(se))

df <- colData(se)
path <- paste("out.csv")
write.csv(df , file = path)
