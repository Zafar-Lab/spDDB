#import anndata as ad
import numpy as np
import scanpy as sc
import pandas as pd
import sys
from scvi.external import RNAStereoscope, SpatialStereoscope

sc_path = sys.argv[1] #Path to single cell data in h5ad format 
st_path = sys.argv[2] #Path to spatial data in h5ad format
annots_key = sys.argv[3] #Path to csv containing cell annotations
output_path = sys.argv[4]
#sc_epochs = sys.argv[4] #Number of epochs for single cell data training
#st_epochs = sys.argv[5] #Number of epochs for spatial data training
sc_data = sc.read_h5ad(sc_path) #Single cell data of dimensions cells x genes
st_data = sc.read_h5ad(st_path) #Spatial data of dimensions spots x genes

sc_data.var_names = sc_data.var_names.astype(str)
st_data.var_names = st_data.var_names.astype(str)
sc_data.var_names_make_unique()
st_data.var_names_make_unique()

sc_data.obs["cell_type"] = sc_data.obs[annots_key]

np.random.seed(1)

## Perform Filtering ##
sc.pp.filter_cells(sc_data, min_genes=1)
sc.pp.filter_genes(sc_data, min_counts = 10) 
non_mito_genes_list = [ name for name in sc_data.var_names if not (name.startswith("MT-") or name.startswith("Mt-"))]

sc_data = sc_data[:, non_mito_genes_list]
sc_data.var_names_make_unique()
sc_data.layers["counts"] = sc_data.X.copy()
sc.pp.normalize_total(sc_data, target_sum = 1e5)
sc.pp.log1p(sc_data)
sc_data.raw = sc_data

sc.pp.highly_variable_genes(
    sc_data,
    n_top_genes = 7000,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    span = 1
)

## Take intersection of genes
intersect = np.intersect1d(st_data.var_names, sc_data.var_names)
sc_data = sc_data[:, intersect].copy()
st_data = st_data[:, intersect].copy()

## Train the scRNA model
RNAStereoscope.setup_anndata(sc_data, layer = "counts", labels_key = "cell_type")
sc_model = RNAStereoscope(sc_data)
sc_model.train(max_epochs = 100)

#Train the spatial model
st_data.layers["counts"] = st_data.X.copy()
SpatialStereoscope.setup_anndata(st_data, layer="counts")
spatial_model = SpatialStereoscope.from_rna_model(st_data, sc_model)
spatial_model.train(max_epochs = 5000)

#Save output
spatial_model.get_proportions().to_csv(output_path)

