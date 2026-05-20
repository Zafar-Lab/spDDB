import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import anndata
import scipy
import seaborn as sns
from matplotlib.pyplot import figure
from sklearn.preprocessing import MinMaxScaler
import torch

np.random.seed(0)
torch.manual_seed(0)


def min_max_scale(vector, min1, max1, min2, max2):
    scaled_vector = ((vector - min1) / (max1 - min1)) * (max2 - min2) + min2
    return np.round(scaled_vector)

def plot_histogram(N_cell2loc):
    
    # plot histogram of no of cells in each spot.
    min1, max1 = np.min(N_cell2loc), np.max(N_cell2loc)
    min2, max2 = 5, 15  # Example target range
    
    scaled_vector = min_max_scale(N_cell2loc, min1, max1, min2, max2)
    
    ##### Plot ################
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (15, 8))
    
    ax1.set_title("N")
    ax1.hist(scaled_vector)
    
    ax2.set_title("N")
    ax2.hist(N_cell2loc)

    return scaled_vector

def generate_expression(ann, sam, original_st_data, decon_matrix, cells_at_spot, ctp_var_names, m_g, spot_barcodes):
    
    decon_matrix = pd.DataFrame(decon_matrix, columns = ctp_var_names)
    
    print (ann.shape)
    print (ann.varm["means_per_cluster_mu_fg"].shape)
    signature_matrix = ann.varm["means_per_cluster_mu_fg"]
    means = np.array(signature_matrix)

    
    ann_X = []
        
    for spot in range(len(decon_matrix)):
    
        col = decon_matrix.columns[0]
        ctp = decon_matrix.iloc[spot][col]
    
        #print (ctp)
        N = cells_at_spot[spot]
    
        N_ct = np.ceil(N*ctp)
    
        if (ctp > 0.0):
            frac_sf = (N*ctp)/N_ct
        else:
            frac_sf = 1.0
    
        gene_exp = (N_ct*signature_matrix["means_per_cluster_mu_fg_" + col]*frac_sf)*m_g
    
        for col in decon_matrix.columns[1:]:
            
            ctp = decon_matrix.iloc[spot][col]
    
            if (ctp > 0.0):
                N_ct = np.ceil(N*ctp)
        
                frac_sf = (N*ctp)/N_ct   
                #print ("ctp", ctp, "frac_sf", frac_sf)
                # Using signature matrix mean representation
                gene_exp += ((N_ct*signature_matrix["means_per_cluster_mu_fg_" + col])*frac_sf)*m_g
    
        
        ann_X += [gene_exp]


    adata = anndata.AnnData(np.array(ann_X))
    
    adata.obs_names = spot_barcodes
    var_names = signature_matrix[signature_matrix.columns[0]].index
    adata.var_names = var_names
    
    for i in original_st_data.uns.keys():
        adata.uns[i] = original_st_data.uns[i]
    
    for i in original_st_data.obs.keys():
        adata.obs[i] = original_st_data.obs[i]
    
    for i in original_st_data.obsm.keys():
        adata.obsm[i] = original_st_data.obsm[i]

    return adata

def compute_distane_matrix(num_spots, tissue_positions_list):
    distance_matrix = np.zeros((num_spots,num_spots))

    for i in range(num_spots):
        
        for j in range(i,num_spots):
            
            x1 = tissue_positions_list[i,0]
            y1 = tissue_positions_list[j,0]
            x2 = tissue_positions_list[i,1]
            y2 = tissue_positions_list[j,1]
            
            distance_matrix[i][j] = (x1 - y1) ** 2 + (x2 - y2) ** 2
            distance_matrix[j][i] = distance_matrix[i][j]
            
    distance_matrix = distance_matrix/np.median(distance_matrix)
    l = 0.1
    #cov_matrix = np.exp(-distance_matrix/l)
    distance_matrix = distance_matrix/np.median(distance_matrix)

    return distance_matrix
    
def sample_UMI(num_spots, distance_matrix, original_st_data, cov_matrix, scale_factor):

    l = 0.1
    sampled_umi = []
    cov_matrices=[]
    
    original_umi_counts = np.array(original_st_data.X.sum(axis=1)).reshape(-1)
     
    for i in range(num_spots):
        for j in range(num_spots):
            
            cov_matrix[i][j] = np.exp(-distance_matrix[i][j]/l)
    
    sampled_umi = np.random.multivariate_normal(original_umi_counts, scale_factor * cov_matrix, check_valid = 'warn')
    
    cov_matrices.append(cov_matrix)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (15,8))
    ax1.set_title("Actual UMI")
    ax1.hist(original_umi_counts)
    
    ax2.set_title("Sampled UMI")
    ax2.hist(sampled_umi)

    factors_change = (sampled_umi/original_umi_counts)
    
    print (np.min(factors_change), np.max(factors_change))
    return sampled_umi, original_umi_counts


    
def scale_anndata(adata, original_st_data, sampled_umi, original_umi_counts):
    sum_simulated_data = np.array(adata.X.sum(axis = 1)).reshape(-1) + 1
    scaling_factors = (sampled_umi/sum_simulated_data)
    
    adata.obs["scaling_factors"] = scaling_factors.reshape(-1)
    print ("scaling ", np.min(scaling_factors), np.max(scaling_factors))
    print (np.any(scaling_factors < 0))
    
    sampled_umi_count = np.array(adata.X.sum(axis = 1)).reshape(-1)
    adata.obs["umi_count_before_scaling"] = sampled_umi_count
    
    # Reshape the scaling factors to a column vector
    scaling_factors = adata.obs["scaling_factors"].values
    
    scaling_factors_column = scaling_factors[:, np.newaxis]
    adata.X = adata.X * scaling_factors_column ## multiply(scaling_factors_column).toarray()
    
    sampled_umi_count = np.array(adata.X.sum(axis = 1)).reshape(-1)
    adata.obs["umi_count_after_scaling"] = sampled_umi_count
    
    original_st_data.obs['umi_counts_actual'] = original_umi_counts
    original_st_data.obs["umi_sampled"] = sampled_umi
    original_st_data.obs["diff"] = (original_umi_counts - sampled_umi)
    adata.obs["diff"] = (adata.obs["umi_count_before_scaling"] - adata.obs["umi_count_after_scaling"])
    
    sc.pl.embedding(original_st_data, basis = "spatial", color = ["umi_counts_actual", "umi_sampled"], cmap='Reds', 
                    save = "original umi.png")
    sc.pl.embedding(original_st_data, basis = "spatial", color = ["diff"], cmap='Spectral', 
                    save = "original umi_spectral.png")
    
    sc.pl.embedding(adata, basis = "spatial", color = ["umi_count_before_scaling", "umi_count_after_scaling"], cmap='Reds', 
                   save = "simulated umi.png")
    sc.pl.embedding(adata, basis = "spatial", color = ["diff"], cmap='Spectral', 
                   save = "simulated umi_spectral.png")

    return adata


















