library("SPOTlight")
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library("anndata")

st_path <- "./simulated_gene_expression_average.h5ad"

sc_path <- "./SC.h5ad"
celltype_key <- 'cell_type'
output_path <- "output.csv"

st<-read_h5ad(st_path)
st$var_names <- make.unique(st$var_names)
st<- as.matrix(st)
st<- t(st)
st <- round(st)
st<- SpatialExperiment(list(counts=st))

sc<-read_h5ad(sc_path)
sc$var_names <- make.unique(sc$var_names)
celltypes <- sc$obs[celltype_key]
sc <- as.matrix(sc)
sc<- t(sc)
sc<- SingleCellExperiment(list(counts=sc))
sc <- logNormCounts(sc)

genes_human <- !grepl(pattern = "^RP[l|s]|MT", x = rownames(sc))
genes_mice <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sc))

dec <- modelGeneVar(sc, subset.row = (genes_human & genes_mice))
hvg <- getTopHVGs(dec, n = 3000)

colLabels(sc) <- celltypes
mgs <- scoreMarkers(sc, subset.row = genes_human & genes_mice, groups = celltypes[,1])

mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  x <- x[x$mean.AUC > 0.8, ]
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)

idx <- split(seq(ncol(sc)), celltypes)
n_cells <- 100
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
})
sc <- sc[, unlist(cs_keep)]

res <- SPOTlight(
  x = sc,
  y = st,
  groups = colLabels(sc)[,1],
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")

final <- res$mat[,sort(colnames(res$mat))]
write.csv(res$mat,output_path)
