
library(STdeconvolve)
library(anndata)
library(scater)

set.seed(101)

#Read SC data
sc_path = "./SC/SN_SC.h5ad"
SC <- read_h5ad(sc_path)
k <- length(unique(SC$obs$cell_type))

#Read ST data
h5adpath <- "./simulated_gene_expression_average.h5ad"
dlpfc <- read_h5ad(h5adpath)
cd <- t(dlpfc$X)

cd <- round(cd)

pos <- dlpfc$obsm$spatial

row.names(pos)<-colnames(cd)
colnames(pos) <- c('x','y')

counts <- cleanCounts(cd)

## feature select for genes
corpus <- restrictCorpus(counts, removeAbove = 1, removeBelow = 0.05)

corpus <- corpus[ ,colSums2(corpus)>0 ]

## choose optimal number of cell-types
ldas <- fitLDA(t(as.matrix(corpus)), Ks = c(k))

## get best model results
optLDA <- optimalModel(models = ldas, opt = "min")

## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta


path <- "output.csv"
write.csv(deconProp , path)


