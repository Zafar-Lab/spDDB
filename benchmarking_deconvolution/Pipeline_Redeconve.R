library(Redeconve)
library("anndata")
library("Matrix")
library("MatrixExtra")
library("dplyr")
set.seed(1)

#memory.limit(size = 400*1024) # in Megabytes (edited) 

args<-commandArgs(trailingOnly = TRUE)
st_path = args[2]
sc_path <-args[1]
celltype_key <- args[3]
output_path <- args[4]

sc<-read_h5ad(sc_path)
sc$var_names <- make.unique(sc$var_names)
celltypes <- sc$obs[celltype_key]
celltypes <- data.frame(barcodes = rownames(celltypes), celltypes = celltypes[,1])
rownames(celltypes) <- celltypes[,1]
sc <- as.matrix(sc)
sc<- t(sc)
sc <- as.csc.matrix(sc)

st<-read_h5ad(st_path)
st$var_names <- make.unique(st$var_names)
st <- as.matrix(st)
st<- t(st)
st <- as.csc.matrix(st)

#Sample 500 cells of each celltype
idx = cell.sampling(ncells = dim(celltypes)[1], celltypes, size = 500, prot=T)
sc.ds = sc[,idx[,1]]
res.ds = deconvoluting(sc.ds,st,genemode="def",hpmode="auto",dopar=T,ncores=16)
final_res <- sc2type(res.ds,celltypes[idx[,1],])
df <- to.proportion(final_res)
df<-t(df)
write.csv(df,output_path)
