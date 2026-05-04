library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(anndata)
library(scran)
options(matrixStats.useNames.NA ="deprecated" )
set.seed(101)

ST <- read_h5ad("ST.h5ad")
q <- length(unique(ST$obs$ground_truth))

x <- ST$obs$array_row
y <- ST$obs$array_col
location <- data.frame('row' = x , 'col' = y)
rownames(location) = colnames(t(ST))

ST <- SingleCellExperiment(assays = list(counts = t(ST$X)) , colData = location , rowData = ST$var)
ST <- logNormCounts(ST)

ST <- spatialPreprocess(ST, platform="Visium")

# Number of clusters
d <- 15  # Number of PCs

## Run BayesSpace clustering
ST <- spatialCluster(ST, q=q, d=d, platform='Visium',
                     nrep=10000, gamma=3, save.chain=TRUE)

#saving results
output <- data.frame("cluster" = ST$spatial.cluster)
rownames(output) <- colnames(ST)
write.csv(output , file ="out.csv")
