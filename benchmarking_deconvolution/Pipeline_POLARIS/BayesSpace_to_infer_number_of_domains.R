
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(anndata)
library(scran)
options(matrixStats.useNames.NA ="deprecated" )
set.seed(101)

start <- Sys.time()

dlpfc <- read_h5ad("./simulated_gene_expression_average.h5ad")


counts <- t(dlpfc$X)
counts <- round(counts)
index <- colSums(counts)>0
counts <- counts[,colSums(counts)>0]

x <- dlpfc$obsm$spatial[index,1]
y <- dlpfc$obsm$spatial[index,2]
location <- data.frame('row' = x , 'col' = y)
rownames(location) = colnames(counts)
                

dlpfc <- SingleCellExperiment(assays = list(counts = counts) , colData = location , rowData = dlpfc$var)

dlpfc <- logNormCounts(dlpfc)

dlpfc <- spatialPreprocess(dlpfc, platform="Visium")

# Number of clusters
d <- 15  # Number of PCs

## Run BayesSpace clustering

qs <- qTune(dlpfc, platform = "Visium", d = d,qs = seq(5,12))
qPlot(qs)

# Number of clusters
q <- 9# based on plot
dlpfc <- spatialCluster(dlpfc, q=q, d=d, platform='Visium',
                        nrep=10000, gamma=3, save.chain=TRUE)


#saving results
final <- data.frame("cluster" = dlpfc$spatial.cluster)
rownames(final) <- colnames(dlpfc)

write.table(final, file ="./POLARIS/simulated_ST_labels.tsv", quote = FALSE , sep = "\t")
  



