library("Matrix")
library("anndata")
library(spacexr)
set.seed(1)

#memory.limit(size = 400*1024) # in Megabytes (edited) 

## Read command line arguments
args<-commandArgs(trailingOnly = TRUE)
sc_path <-args[1]
st_path <- args[2]
celltype_key <- args[3]
output_path <- args[4]

st<-read_h5ad(st_path)
st$var_names <- make.unique(st$var_names)

## Read spatial coordinates and create ST object
new_x <- st$obsm$spatial[,2]
new_y <- st$obsm$spatial[,1]
#if ("array_row" %in% colnames(st$obs)) {
#	new_x <- st$obs['array_row']
#	new_y <- st$obs['array_col']
#} else if("new_x" %in% colnames(st$obs)){
#	new_x <- st$obs['new_x']
#	new_y <- st$obs['new_y']
#} else {
#  new_x <- st$obsm$spatial[,2]
#  new_y <- st$obsm$spatial[,1]
#}
coords_df <- data.frame(cbind(new_x,new_y))
rownames(coords_df) <- rownames(st)
colnames(coords_df) <- cbind('x','y')
st<- as.matrix(st)
st<- t(st)
#st <- round(st)
spatial <- SpatialRNA(coords_df,st)

## Read sc data and create SC object
sc<-read_h5ad(sc_path)
sc <- sc[rowSums(sc$X)!=0] # Ensure no empty cells
sc$var_names <- make.unique(sc$var_names)
celltypes <- sc$obs[celltype_key]
celltypes[,1] <- gsub("/", "-", celltypes[,1]) #'/' not allowed in celltype names

## Remove celltypes with counts less than 25 and convert celltypes to correct format
bool <- celltypes[,1] %in% names(table(celltypes)[table(celltypes)>=25])
removed_ct <- names(table(celltypes)[table(celltypes)<25])
sc <- sc[bool]
celltypes <- sc$obs[celltype_key]
celltypes[,1] <- gsub("/", "-", celltypes[,1])
celltypes_factor <- as.factor(celltypes[,celltype_key])
names(celltypes_factor) <- rownames(celltypes)

sc <- as.matrix(sc)
sc<- t(sc)
reference <- Reference(sc, celltypes_factor)

## Run RCTD deconvolution
myRCTD <- create.RCTD(spatial, reference, max_cores = 4, UMI_min=100)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

##Get output in correct format
weights <- myRCTD@results$weights
norm_weights <- normalize_weights(weights)
norm_weights <- as.data.frame(as.matrix(norm_weights))

## Add 0s for those celltypes that were removed
if(length(removed_ct)>0){
  for (i in 1:length(removed_ct)){
    norm_weights[removed_ct[i]] = rep(0,dim(norm_weights)[1])
  }
}
norm_weights<-norm_weights[,order(colnames(norm_weights))]
write.csv(norm_weights, output_path)
