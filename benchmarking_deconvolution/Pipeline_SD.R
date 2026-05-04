setwd("./SD/SD2")
source("SD2_utiles.R")
library(anndata)
library(scater)


###helper functions
generate_spots = function (se_obj, clust_vr,lower = 2, upper = 10, n = 1000, verbose = TRUE) 
{
  if (is(se_obj) != "Seurat") 
    stop("ERROR: se_obj must be a Seurat object!")
  if (!is.character(clust_vr)) 
    stop("ERROR: clust_vr must be a character string!")
  if (!is.numeric(n)) 
    stop("ERROR: n must be an integer!")
  if (!is.logical(verbose)) 
    stop("ERROR: verbose must be a logical object!")
  suppressMessages(require(DropletUtils))
  suppressMessages(require(purrr))
  suppressMessages(require(dplyr))
  suppressMessages(require(tidyr))
  se_obj@meta.data[, clust_vr] <- gsub(pattern = "[[:punct:]]|[[:blank:]]", 
                                       ".", x = se_obj@meta.data[, clust_vr], perl = TRUE)
  print("Generating synthetic test spots...")
  start_gen <- Sys.time()
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  count_mtrx <- as.matrix(se_obj@assays$RNA@counts)
  ds_spots <- lapply(seq_len(n), function(i) {
    cell_pool <- sample(colnames(count_mtrx), sample(x = lower:upper, 
                                                     size = 1))
    pos <- which(colnames(count_mtrx) %in% cell_pool)
    
    #### EDITED HERE, ADDED drop = FALSE
    tmp_ds <- se_obj@meta.data[pos, , drop= FALSE ] %>% mutate(weight = 1)
    name_simp <- paste("spot_", i, sep = "")
    #### EDITED HERE, ADDED SUBCLASS TO SELECT()
    spot_ds <- tmp_ds %>% dplyr::select(all_of(clust_vr), subclass ,
                                        weight ) %>% dplyr::group_by(!!sym(clust_vr)) %>% 
      dplyr::summarise(sum_weights = sum(weight)) %>% 
      dplyr::ungroup() %>% tidyr::pivot_wider(names_from = all_of(clust_vr), 
                                              values_from = sum_weights) %>% dplyr::mutate(name = name_simp)
    syn_spot <- rowSums(as.matrix(count_mtrx[, cell_pool]))
    sum(syn_spot)
    names_genes <- names(syn_spot)
    if (sum(syn_spot) > 25000) {
      syn_spot_sparse <- DropletUtils::downsampleMatrix(Matrix::Matrix(syn_spot, 
                                                                       sparse = T), prop = 20000/sum(syn_spot))
    }
    else {
      syn_spot_sparse <- Matrix::Matrix(syn_spot, sparse = T)
    }
    rownames(syn_spot_sparse) <- names_genes
    colnames(syn_spot_sparse) <- name_simp
    setTxtProgressBar(pb, i)
    return(list(syn_spot_sparse, spot_ds))
  })
  ds_syn_spots <- purrr::map(ds_spots, 1) %>% base::Reduce(function(m1, 
                                                                    m2) cbind(unlist(m1), unlist(m2)), .)
  ds_spots_metadata <- purrr::map(ds_spots, 2) %>% dplyr::bind_rows() %>% 
    data.frame()
  ds_spots_metadata[is.na(ds_spots_metadata)] <- 0
  lev_mod <- gsub("[\\+|\\ |\\/]", ".", unique(se_obj@meta.data[,clust_vr]))
  colnames(ds_spots_metadata) <- gsub("[\\+|\\ |\\/]", ".", 
                                      colnames(ds_spots_metadata))
  print(sum(lev_mod %in% colnames(ds_spots_metadata)))
  print(length(unique(se_obj@meta.data[,clust_vr])) + 1)
  if (sum(lev_mod %in% colnames(ds_spots_metadata)) == (length(unique(se_obj@meta.data[,clust_vr])) + 1)) {
    ds_spots_metadata <- ds_spots_metadata[, lev_mod]
  }
  else {
    missing_cols <- lev_mod[which(!lev_mod %in% colnames(ds_spots_metadata))]
    ds_spots_metadata[missing_cols] <- 0
    ds_spots_metadata <- ds_spots_metadata[, lev_mod]
  }
  close(pb)
  print(sprintf("Generation of %s test spots took %s mins", 
                n, round(difftime(Sys.time(), start_gen, units = "mins"), 
                         2)))
  print("output consists of a list with two dataframes, this first one has the weighted count matrix and the second has the metadata for each spot")
  return(list(topic_profiles = ds_syn_spots, cell_composition = ds_spots_metadata))
}


gen_pseudo_ST <- function(st_count,st_label,spot_num = 1000,
                          HVG_num = 200,scale_num = 10000,lower_cellnum,upper_cellnum){
  
  feature1 = rownames(extract_dropout_genes(st_count[[1]]))
  print(length(feature1))
  feature2 = select_feature(st_count[[1]],st_label[[1]],nf = HVG_num)
  print(length(feature2))
  print(dim(st_count[[2]]))
  sel.features = union(feature1,feature2)
  print(length(sel.features))
  
  #edited here to resolve NA being introduced as feature for imaging data
  sel.features <- na.omit(sel.features)
  
  st_count_new <- list(st_count[[1]][sel.features,],st_count[[2]][sel.features,])
  print(sprintf('dropout_features:%s',length(feature1)))
  print(sprintf('HVG_features:%s',length(feature2)))
  
  colnames(st_label[[1]]) <- 'subclass'
  #### EDITED HERE, CREATED SEURAT OBJECT INDIRECTY THROUGH SCE TO GET APPROPRIATE FORM
  # tem.t1 <- Seurat::CreateSeuratObject(counts = st_count_new[[1]],meta.data=st_label[[1]]);
  
  temp <- SingleCellExperiment(assays = list(counts = st_count_new[[1]]))
  temp <- logNormCounts(temp)
  temp2 <- as.Seurat(temp , counts = "counts")
  temp2 <- RenameAssays(temp2 , "originalexp" , "RNA")
  temp2@meta.data <- st_label[[1]]
  tem.t1 <- temp2
  #' convert scRNA-seq data to pseudo-spatial data                                                                                                                           
  test.spot.ls1<-generate_spots(se_obj=tem.t1,clust_vr='subclass',n = spot_num,lower = lower_cellnum,upper = upper_cellnum);
  test.spot.counts1 <- as.matrix(test.spot.ls1[[1]])
  colnames(test.spot.counts1)<-paste("mixt",1:ncol(test.spot.counts1),sep="_");
  test.spot.metadata1 <- test.spot.ls1[[2]]
  
  st_counts <- list(test.spot.counts1,st_count_new[[2]])
  st_labels <- list(test.spot.metadata1/rowSums(test.spot.metadata1))
  st_norm <- normalize_data(st_counts,scale_num = scale_num)[[1]]
  st_scale <- scale_data(st_counts,st_norm,sel.features)
  return (list(st_counts,st_labels,st_norm,st_scale,sel.features))
}
set.seed(101)

st_anndata<- read_h5ad("./simulated_gene_expression_average.h5ad")#//genes*spots
sc_anndata<- read_h5ad("./SC/SN_SC.h5ad")##genes*cells
sc_matrix<- t(sc_anndata$X)
st_matrix <- t(st_anndata$X)

cell_type <- sc_anndata$obs$cell_type

spatial_location_anndata<- data.frame(x=st_anndata$obsm$spatial[,1], y=st_anndata$obsm$spatial[,2])

row.names(spatial_location_anndata)<-colnames(st_matrix)

sc_matrix <- as.matrix(sc_matrix)
st_matrix <- as.matrix(st_matrix)

index <- colSums(sc_matrix) > 0
sc_matrix <- sc_matrix[,index]
cell_type = cell_type[index]

SD2(sc_matrix,
    st_matrix,
    cell_type,
    ST_location = spatial_location_anndata,
    spot_num = 300, 
    lower_cellnum = 10,
    upper_cellnum = 20)

##########################################################
## run ./SD2/train.py before running further lines of code

setwd("./SD/SD2")
output <- read.csv("./SD2_Result/predict_output.csv")
rownames(output) <- rownames(st_anndata)
output <- output[,sort(colnames(output))]

path <- "output.csv"
write.csv(output , path)

