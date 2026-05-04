#If version conflicts between anndata and giotto python environments:
# .rs.restartR()
# 
# reticulate::use_python("/usr/bin/python3.8", required = TRUE)

library(Giotto)
library(DescTools)
library(anndata)
st_path <- "./ST.h5ad"
ST<- read_h5ad(st_path)
ST$var_names_make_unique()
st_adata <- t(as.matrix(ST$X))

x<- ST$obs$array_row
y<- ST$obs$array_col
spatial_location_anndata<- data.frame(x=x, y=y)

row.names(spatial_location_anndata)<-colnames(st_adata)

#Refresh R to restart python env if conflicts
# .rs.restartR()
# 
# reticulate::use_condaenv(required = T,  '/home/user/.local/share/r-miniconda/envs/giotto_env/bin/python')
library(Giotto)
library(DescTools)
library(anndata)

my_giotto_object = createGiottoObject(raw_exprs = st_adata,
                                      spatial_locs = spatial_location_anndata)

# processing
visium_brain <- filterGiotto(gobject = my_giotto_object,
                             expression_threshold = 1,
                             # gene_det_in_min_cells = 50,
                             # min_det_genes_per_cell = 500,
                             expression_values = c('raw'),
                             verbose = T)
visium_brain <- normalizeGiotto(gobject = visium_brain, scalefactor = 6000, verbose = T)

# HVG step removed for certain simulated visium datasets:
visium_brain <- calculateHVG(gobject = visium_brain)

visium_brain <- runPCA(gobject = visium_brain, 
                       scale_unit = T, 
                       center=T, 
                       method="factominer")

#changed to 1:10 if error
visium_brain <- runUMAP(visium_brain, dimensions_to_use = 1:20)
visium_brain <- runtSNE(visium_brain, dimensions_to_use = 1:20)

visium_brain <- createNearestNetwork(gobject=visium_brain, 
                                     
)

#resolution changed as per number of clusters:
my_giotto_object = doLeidenCluster(visium_brain , python = "/home/user/.local/share/r-miniconda/envs/giotto_env/bin/python", resolution = 1.78)
final <- as.data.frame(my_giotto_object@cell_metadata)
length(unique(final$leiden_clus))

path <-"out.csv"
write.csv(final,path)

