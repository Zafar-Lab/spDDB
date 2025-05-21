library(reticulate)
library(anndata)
library(zellkonverter)
library(SummarizedExperiment)
library(Matrix)

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      crop = NA)

## ----warning = FALSE----------------------------------------------------------
suppressPackageStartupMessages({
  library(countsimQC)
  library(DESeq2)
})


convert_ann2dds_original <- function(orig_path, sim_path) {
    ann_orig <- read_h5ad(orig_path)
    ann_sim <- read_h5ad(sim_path)

    genes_orig <- ann_orig$var_names
    genes_sim <- ann_sim$var_names

    common_genes <- intersect(genes_orig, genes_sim)
    obs_orig <- ann_orig$obs_names
    obs_sim <- ann_sim$obs_names
    cells_to_filter <- intersect(obs_orig, obs_sim)
    
    ann_filtered <- ann_orig[ , colnames(ann_orig) %in% common_genes]
    ann_filtered <- ann_filtered[ann_filtered$obs_names %in% cells_to_filter, ]
    
    ann <- ann_filtered

    print (ann$shape)
    print (ann_sim$shape)
    #print (common_genes$shape)
    
    #print (ann$layers["raw"])
    counts <- ann$X #as.matrix(ann$X)
    rownames(counts) <- make.names(seq_len(nrow(counts)), unique = TRUE)
    colnames(counts) <- make.names(seq_len(ncol(counts)), unique = TRUE)
    attributes(rownames(counts)) <- NULL
    attributes(colnames(counts)) <- NULL
    
    
    if (is(counts, "sparseMatrix")) {
        message("The matrix is sparse")
        counts <- as.matrix(counts)
        
    } else {
        message("The matrix is already dense")
    }
    
    counts <- t(counts)
    metadata <- as.data.frame(ann$obs)           
    rownames(metadata) <- colnames(counts)
    
    #print (all(rownames(metadata) == colnames(counts))
    
    dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ 1
    )     
}

convert_ann2dds <- function(path) {
    ann <- read_h5ad(path)

    #print (ann$layers["raw"])
    counts <- ann$X #as.matrix(ann$X)
    rownames(counts) <- make.names(seq_len(nrow(counts)), unique = TRUE)
    colnames(counts) <- make.names(seq_len(ncol(counts)), unique = TRUE)
    attributes(rownames(counts)) <- NULL
    attributes(colnames(counts)) <- NULL
    
    
    if (is(counts, "sparseMatrix")) {
        message("The matrix is sparse")
        counts <- as.matrix(counts)
        
    } else {
        message("The matrix is already dense")
    }
    
    counts <- t(counts)
    metadata <- as.data.frame(ann$obs)           
    rownames(metadata) <- colnames(counts)
    
    #print (all(rownames(metadata) == colnames(counts))
    
    dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ 1
    )
}

countsim_plots <- function(ddsList, name_to_save){
    tempDir <- getwd()
    countsimQCReport(ddsList = ddsList, outputFile = paste0(name_to_save, ".html"),
                 outputDir = tempDir, outputFormat = "html_document", 
                 showCode = FALSE, forceOverwrite = TRUE,
                 savePlots = TRUE, description = "This is my test report.", 
                 maxNForCorr = 25, maxNForDisp = Inf, 
                 calculateStatistics = TRUE, subsampleSize = 25,
                 kfrac = 0.01, kmin = 5, 
                 permutationPvalues = FALSE, nPermutations = NULL)
        
}

individual_plots <- function(name_to_save) {
    tempDir <- getwd()
    #file.path
    ggplots <- readRDS(paste0(tempDir, paste0("/", paste0(name_to_save, "_ggplots.rds"))))
    
    if (!dir.exists(file.path(tempDir, paste0("figures_", name_to_save)))) {
        dir.create(file.path(tempDir, paste0("figures_", name_to_save)))
    }
    generateIndividualPlots(ggplots, device = "pdf", nDatasets = 2, 
                            outputDir = file.path(tempDir, paste0("figures_", name_to_save)))    
}

