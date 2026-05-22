import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from diptest import diptest

# Getting the anndata between coordinates (x1, y1) and (x2, y2)

# This is giving a good view of how the axis is shifting for the binning purposes.

def plot_lines(coordinates, x1, y1, x2, y2):
    
    plt.scatter(coordinates[:, 0], coordinates[:, 1], s=1)
    plt.axvline(x = x1, color = 'red')
    plt.axvline(x = x2, color = 'red')
    plt.axhline(y = y1, color = 'red')
    plt.axhline(y = y2, color = 'red')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Spatial coordinates with bounding box')
    plt.show()
    
def getspots(adata, x1, y1, x2, y2, coord1, coord2):

    coordinates = adata.obsm['spatial']
    #plot_lines(coordinates, x1, y1, x2, y2)
    
    mask = (
        (coordinates[:, coord1] >= x1) & (coordinates[:, coord1] <= x2) &
        (coordinates[:, coord2] >= y1) & (coordinates[:, coord2] <= y2)
    )

    #if (np.sum(mask) != 0):
    #    print(f"Number of spots in mask: {np.sum(mask)}")
        
    adata_subset = adata[mask].copy()
    #print ("(",x1, ",", y1, ") to (", x2, ",", y2, ")", adata_subset)
    
    return adata_subset

def binning (unbinned_st_path, spot_diameter, celltype_col, dataset, coord1, coord2):

    
    st_adata = sc.read_h5ad(unbinned_st_path)
    
    if (dataset in ["merfish_brain", "merfish_breast_cancer"]):
        print (dataset)
        st_adata.obsm["spatial"] = st_adata.obsm["spatial"].to_numpy()

    #print (st_adata.obsm["spatial"])
    cell_types = st_adata.obs[celltype_col].unique().to_list()
    
    min_x = st_adata.obsm["spatial"][:, coord1].min()
    max_x = st_adata.obsm["spatial"][:, coord1].max()
    min_y = st_adata.obsm["spatial"][:, coord2].min()
    max_y = st_adata.obsm["spatial"][:, coord2].max()

    
    
    ann_list = []
    
    itr_y = min_y
    while (itr_y <= max_y - spot_diameter):
        
        itr_x = min_x
        
        while (itr_x <= max_x - spot_diameter):

            cell_type_dict = {cell_type: 0 for cell_type in cell_types}
            ann_subset = getspots(st_adata, itr_x, itr_y, itr_x + spot_diameter, itr_y + spot_diameter, coord1, coord2)
    
            if (ann_subset.shape[0] != 0 ):
                x_row = np.array(ann_subset.X.sum(axis = 0).reshape(1, -1)) # row wise sum of all the selected spots.
                #print (x_row)
                obs_row = ann_subset.obs[celltype_col].value_counts().to_dict()
                spatial_coords = [ann_subset.obsm["spatial"][:, 0].mean(), ann_subset.obsm["spatial"][:, 1].mean()]

                for key in obs_row:
                    cell_type_dict[key] = obs_row[key]
    
                no_of_cells = pd.DataFrame(cell_type_dict, index=[0]).astype(int)
                #print (no_of_cells)
                binned_ann = anndata.AnnData(X = x_row, obs = no_of_cells)
                binned_ann.obsm["spatial"] = np.array([spatial_coords])
                ann_list += [binned_ann]
                
                #break
            itr_x += spot_diameter
        
        #break
        itr_y += spot_diameter

    combined_adata = sc.concat(ann_list, join='outer', label='dataset')
    return combined_adata

def diptest(df):
    
    for column in df.columns:
        sns.kdeplot(df[column], label=column, bw_adjust=0.5)
    
    plt.title('KDE Plots')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.legend(title ='Columns')
    plt.show()
