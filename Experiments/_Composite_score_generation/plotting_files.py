import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import warnings
warnings.filterwarnings("ignore")

lst_columns = ['id', 
                # Composite scores....
                'overall_composite_score', 'overall_composite_score_label', 'comp_acc', 'comp_acc_label', 
                'comp_bivariate', 'comp_bivariate_label', 'composite_cellcharter_score', 'composite_cellcharter_score_label',
                'composite_rare_score', 'composite_rare_score_label',
                # Spatial metrics...
                'morans_r', 'morans_r_label', 'lee_stat', 'lee_stat_label',  
                'spatial_pearson', 'spatial_pearson_label', 'ssim', 'ssim_label', 'geary_c', 'geary_c_label',         
                # Non-spatial metrics...
                'cosine_sim', 'cosine_sim_label', 'pearson', 'pearson_label',  'JSD', 'JSD_label', 'RMSE', 'RMSE_label', 
                # CellCharter metrics
                'comp_curl_top_spatial', 'comp_curl_top_spatial_label',
                'comp_elongation_top_spatial', 'comp_elongation_top_spatial_label', 'comp_linearity_top_spatial', 
                'comp_linearity_top_spatial_label',
                'comp_elongation_bottom_spatial', 'comp_elongation_bottom_spatial_label', 
                'comp_linearity_bottom_spatial', 'comp_linearity_bottom_spatial_label',                                    
                #Rare cell type metrics
                'comp_rare_def1_spatial', 'comp_rare_def1_spatial_label', 'comp_rare_def2_spatial', 'comp_rare_def2_spatial_label', 
               "AUPR", "AUPR_label"
          ]

lst_group = ["Method"] + \
                     ["Composite Score" for i in range(10)] + \
                    ["Spatial metrics" for i in range(10)] + \
                    ["Non-Spatial metrics" for i in range(8)] + \
                    ["Cell Charter Metrics" for i in range(10)] + \
                    ["Rare celltype Metrics" for i in range(6)]

def create_columns_group():    

    
    columns_group = pd.DataFrame(columns = ["Experiment", "Category", "group", "palette"])
    columns_group["group"] = lst_group
    columns_group["Experiment"] = lst_group
    columns_group["Category"] = lst_group #["Method"] + ["celltype level metric" for i in range(32)]
    
    columns_group["palette"] = ["overall"] + \
                                ["comp_score" for i in range(10)] + \
                                ["score" for i in range(10)] + \
                                ["accuracy" for i in range(8)] + \
                                ["cell_charter" for i in range(10)] + \
                                 ["rare" for i in range(6)]
    
    return columns_group

def create_columns_info(methods_all):
    
    columns = pd.DataFrame(columns = ["group", "id", "name", "geom"])

    lst_names = ["Method", 
    
                 "Overall Composite Score", pd.NA,
                 "Composite Non-Spatial Score", pd.NA, "Composite Spatial Score", pd.NA,
                 "Composite CellCharter Score", pd.NA, "Composite Rare Score", pd.NA,
                 
                 "Moran's I", pd.NA, "Lee's L statistics", pd.NA,  "Spatial Pearson Correlation", 
                 pd.NA, "SSIM", pd.NA,  "Geary' C", pd.NA, "Cosine Similarity", pd.NA, "Pearson Correlation", pd.NA,
                "JS", pd.NA, "RMSE", pd.NA, 
                 
                 
                "Composite curl score-top%", pd.NA, "Composite elongation score-top%", pd.NA, "Composite linearity score-top%", pd.NA,
                "Composite elongation score-bottom%", pd.NA, "Composite linearity score-bottom%", pd.NA, 
                 
                 "Composite rare score-def1", pd.NA, "Composite rare score-def2", pd.NA,
                 "AUPR", pd.NA
                ]
    
    lst_type = ["text"] + ["bar", "text"]*5 + ["circle", "text"]*9 + ["circle", "text"]*8
    
    columns["id"] = lst_columns
    columns["group"] = lst_group
    
    columns["name"] = lst_names
    columns["geom"] = lst_type
    columns["size"] = [4] + [2, 3]*22
    
    return columns


def create_rows_group_and_info(methods_all):
    row_groups = pd.DataFrame(columns = ["group", "Group"])

    lst = ["Linear-model", "Bayesian-model", "Bayesian-model", "Deep-learning-model", "Deep-learning-model", 
                    "Deep-learning-model", "Deep-learning-model", "Linear-model", "Linear-model", "Graph-based-model", "Anchor-based-model",
                     "Bayesian-model",  "Bayesian-model", "NMI-based-model", "w/o scRNA-seq Ref", "NMI-based-model", "Linear-model",
                    "Bayesian-model", "Deep-learning-model", "Bayesian-model", "Deep-learning-model"]
    
    #lst = ["." for l in lst]    
    row_groups["group"] = lst
    row_groups["Group"] = lst       
    rows = pd.DataFrame(columns = ["group", "id"])    
    rows["id"] = methods_all #grouped_df["id"].values
    rows["size"] = [2 for i in range(rows.shape[0])]       
    rows["group"] = lst    
    return row_groups, rows



