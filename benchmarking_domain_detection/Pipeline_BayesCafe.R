setwd("./BayesCafe-main")
source("R/BayesCafe.R")

library(anndata)
library(scater)
set.seed(101)

ST <- read_h5ad("./ST.h5ad")
out.path <- "./output_BayesCafe.csv"
ST$obs_names <- make.unique(ST$obs_names)
ST$var_names <- make.unique(ST$var_names)

count_ST <- (ST$X)

x<- ST$obsm$spatial[,1]
y<- ST$obsm$spatial[,2]
pos <- data.frame('x' = x , 'y' = y)
rownames(pos) <- rownames(ST)

result <- dataPreprocess(
  count = count_ST, 
  loc = pos,
  cutoff_sample = 100,
  cutoff_feature = 0.1,
  cutoff_max = 0,
  size.factor = "tss",
  platform = "Visium",
  findHVG = FALSE,
  n.HVGs=2000)

count <- result$count
loc <- result$loc
s <- result$s
P <- result$P

res <- bayes_cafe(
  count = count, 
  K = 12, 
  s = s, 
  P = P,
  iter = 1200,
  burn = 1000
)

write.csv(res$cluster_result , out.path)