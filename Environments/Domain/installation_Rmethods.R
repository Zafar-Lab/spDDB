install.packages("BioManager")
BiocManager::install("SPOTlight")
BiocManager::install("SingleCellExperiment")
BiocManager::install("SpatialExperiment")
BiocManager::install("scater")
BiocManager::install("scran")
BiocManager::install("SpatialDecon")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

install.packages("devtools")
library(devtools)
install_github("https://github.com/MarcElosua/SPOTlight")
library(spacexr) #RCTD

library(ggplot2)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)

# install.packages("devtools")


devtools::install_github("timothyhyndman/deconvolve")
library(deconvolve)

install.packages('Seurat')
library(Seurat)
#remotes::install_version("Seurat", version = "3.X.X")
library(SpatialDecon)

remotes::install_github("r-spatial/sf")
remotes::install_github("rspatial/raster")
remotes::install_github("mtennekes/tmaptools") # required for dev version of tmap
remotes::install_github("mtennekes/tmap")
library(sf)
library(raster)
library(tmap)

