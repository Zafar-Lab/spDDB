# if(!require(devtools))
#   install.packages(devtools)
# devtools::install_github("zhengli09/BASS")

library(BASS)
library(scater)
library(anndata)
set.seed(101)


for (i in c(151507:151510,151669:151676)){
  if (i %in% c(151507:151510,151673:151676)){
    R <- 7
  } else {
    R <- 5
  }
  i <- as.character(i)
  ST <- read_h5ad(paste0("/data/Ajita/Spatial/Datasets/Spatial_Deconvolution/DLPFC/ST/",i,"/ST.h5ad"))
  C <- 20
  # R <- 7
  
  ST$obs_names <- make.unique(ST$obs_names)
  ST$var_names <- make.unique(ST$var_names)
  
  cnts <- list(t(ST$X))
  
  x <- ST$obs$array_row
  y <- ST$obs$array_col
  xy <- data.frame('row' = x , 'col' = y)
  # xy <- ST$obsm$spatial
  rownames(xy) <- rownames(ST)
  xy <- list(xy)
  
  
  
  BASS <- createBASSObject(cnts, xy, C , R , beta_method = "SW")
  
  BASS <- BASS.preprocess(BASS)
  
  BASS <- BASS.run(BASS)
  
  BASS <- BASS.postprocess(BASS)
  
  final <- data.frame("cluster" = BASS@results$z[[1]])
  rownames(final) <- ST$obs_names
  
  path <- paste0("/home/keerthana/Documents/Spatial/spatial_domain/final_outputs/DLPFC/",i,"/output_",i,"_BASS.csv")
  write.csv(final , path)
}
