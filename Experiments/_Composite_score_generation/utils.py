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
from plotting_files import *
import re

replacements = {
    "Mouse_brain_ST8059048" : {"data_name" : "Mouse Brain - ST48", "cc_data_name" : "Mouse_brain_ST48"},
    "Mouse_brain_ST8059052" : {"data_name" : "Mouse Brain - ST52", "cc_data_name" : "Mouse_brain_ST52"},
    "Breast_Atlas_ST_1" : {"data_name" : "Breast_Atlas_ST1", "cc_data_name" : "Breast_Atlas_ST1"},
    "Breast_Atlas_ST_8" : {"data_name" : "Breast_Atlas_ST8", "cc_data_name" : "Breast_Atlas_ST8"},
    "Liver_Atlas_4" : {"cc_data_name" : "Liver_Atlas_Mouse_4", "data_name" : "Liver Atlas Mouse - 4"},
    "Liver_Atlas_12" : {"cc_data_name" : "Liver_Atlas_Mouse_12", "data_name" : "Liver Atlas Mouse - 12"},
    "Liver_Atlas_18" : {"cc_data_name" : "Liver_Atlas_Human", "data_name" : "Liver Atlas Human - 18"},
    "Kidney_Atlas_200903_1" : {"data_name" : "Kidney Atlas - 1", "cc_data_name" : "Kidney_Atlas_Human_1"},
    "Kidney_Atlas_200903_6" : {"data_name" : "Kidney Atlas - 6", "cc_data_name" : "Kidney_Atlas_Human_6"}
}

def filter_cell_types(cc_path, sheet_path, data_name):

    data_name  = str(data_name)
    df_C = pd.read_excel(open(cc_path, 'rb'), sheet_name = "C_or")
    df_E = pd.read_excel(open(cc_path, 'rb'), sheet_name = "E_or")
    df_L = pd.read_excel(open(cc_path, 'rb'), sheet_name = "L_or")

    # Finding curl, linearity and elongation metrics

    cc_data_name = data_name

    if ("DLPFC" in data_name):
        # 151508 is not in str format.
        cc_data_name = int(data_name.split("_")[1])

    # Condition only for cancer datasets
    elif (data_name not in df_C["Dataset"].unique()) and ("Cancer" in data_name):
        print ("replacing space with _")
        cc_data_name = data_name.replace(" ", "_") # For Kidney and Breast Cancer
        cc_data_name = cc_data_name.replace("HCC-", "") # For Liver Cancer
        data_name = cc_data_name

    if (data_name in replacements.keys()):
        print ("key", data_name)
        cc_data_name = replacements[data_name]["cc_data_name"]
        data_name = replacements[data_name]["data_name"]
        
    print (data_name, "and", cc_data_name)
    
    df_C = df_C[df_C["Dataset"] == cc_data_name]
    df_C_top = df_C[df_C["type"] == "top"]
    c_top_ct = [s[:-2] for s in df_C_top["cell_type"].values]
    
    df_E = df_E[df_E["Dataset"] == cc_data_name]
    df_E_top = df_E[df_E["type"] == "top"]
    df_E_bottom =  df_E[df_E["type"] == "bottom"]
    e_top_ct = [s[:-2] for s in df_E_top["cell_type"].values]
    e_bottom_ct = [s[:-2] for s in df_E_top["cell_type"].values]
    
    df_L = df_L[df_L["Dataset"] == cc_data_name]
    df_L_top = df_L[df_L["type"] == "top"]
    df_L_bottom =  df_L[df_L["type"] == "bottom"]
    l_top_ct = [s[:-2] for s in df_L_top["cell_type"].values]
    l_bottom_ct = [s[:-2] for s in df_L_bottom["cell_type"].values]
    
    ### Finding rare cell types
    df_rare2 = pd.read_excel(open(sheet_path, 'rb'), sheet_name = "Rare-Def2")
    df_rare2["cell_types"] = df_rare2["Celltypes-Def 2"].str.split("_count").str[0].str.rstrip().str.replace('-', '_').values
    rare2_cts = df_rare2[df_rare2["Datasets"].isin([data_name])]["cell_types"].values.tolist()
    #print (df_rare2["Datasets"].unique())
    
    df_rare1 = pd.read_excel(open(sheet_path, 'rb'), sheet_name = "Rare-Def1")
    df_rare1["cell_types"] = df_rare1["Celltypes-Def 1"].str.split(".0").str[0].str.rstrip().str.replace('-', '_').values
    rare1_cts = df_rare1[df_rare1["Datasets"].isin([data_name])]["cell_types"].values.tolist()
    #print (df_rare1["Datasets"].unique())
    
    symbols_to_replace = ["-", " ", "/", ".", "+", "(", ")", ",", "&"]

    #symbols_to_replace = ['-', '\\', '.']

    pattern = '[' + re.escape(''.join(symbols_to_replace)) + ']'
    
    # Replace the symbols in each string
    rare1_cts = [x for x in rare1_cts if x is not np.nan]
    rare2_cts = [x for x in rare2_cts if x is not np.nan]
    
    if (len(rare1_cts) > 0):
        rare1_cts = [re.sub(pattern, '_', s) for s in rare1_cts]
    else:
        rare1_cts = []

    if (len(rare2_cts) > 0):
        rare2_cts = [re.sub(pattern, '_', s) for s in rare2_cts]
    else:
        rare2_cts = []
        
    l_top_ct =  [re.sub(pattern, '_', s) for s in l_top_ct]
    l_bottom_ct = [re.sub(pattern, '_', s) for s in l_bottom_ct]
    c_top_ct = [re.sub(pattern, '_', s) for s in c_top_ct]
    e_top_ct = [re.sub(pattern, '_', s) for s in e_top_ct]
    e_bottom_ct = [re.sub(pattern, '_', s) for s in e_bottom_ct]
    
    print ("linearity_top_ct -->", l_top_ct, "\n", \
       "linearity_bottom_ct -->", l_bottom_ct, "\n", \
       "df_C_top --->", c_top_ct, "\n", \
       "elongation_top_ct", e_top_ct, "\n", 
       "elongation_bottom_ct", e_bottom_ct, "\n",
      "rare2_cts", rare2_cts, "\n",
      "rare1_cts", rare1_cts, "\n")

    return l_top_ct, l_bottom_ct, c_top_ct, e_top_ct, e_bottom_ct, rare2_cts, rare1_cts

def plot_minmax(df, methods, all_metrics):
    for mt in all_metrics:
        df_ = df[df["Metric"] == mt]
    
        minimum = np.min(np.min(df_[methods]))
        maximum = np.max(np.max(df_[methods]))
        print (mt, "-->", minimum, maximum)
        
def scaing_of_dataframe(df, methods, all_metrics):
    
    #plot_minmax(df, methods, all_metrics)
    # Below steps were needed for scaling per metric as dataframe has all the metrics
    df_list = []    
    for m in df["Metric"].unique():
        
        df_mt = df[df["Metric"] == m]
        cell_types = df_mt["cell_type"].values
        df_mt = df_mt[methods]
        #df_mt = df_mt[df_mt.columns[4:]]
    
        min_val = df_mt.min().min()
        max_val = df_mt.max().max()
        
        df_scaled = (df_mt - min_val) / (max_val - min_val)
        df_scaled["Metric"] = m
        df_scaled["cell_type"] = cell_types
    
        df_list += [df_scaled]
    
    merged_df = pd.concat(df_list, axis = 0, ignore_index = True)
    #plot_minmax(merged_df, methods, all_metrics)
    return merged_df

# df.mean by default excludes string columns such as celltypes and methods.
def composite_score_rare(grouped_df, df_bivariate, rare1_cts, rare2_cts, methods_all):

    print ("CELL TYPES", df_bivariate["cell_type"].unique(), "COLUMNS are", df_bivariate.columns)
    df_rare1 = df_bivariate[df_bivariate["cell_type"].isin(rare1_cts)].reset_index(drop = True)

    rare1_mean = df_rare1[methods_all]
    
    print ("Main dataframe", grouped_df.shape, rare1_mean.shape, "==", rare1_mean.mean().values, len(rare1_mean.mean().values,))

    grouped_df["comp_rare_def1_spatial"] = rare1_mean.mean(axis=0).values
    
    if len(rare2_cts) > 0:

        df_rare2 = df_bivariate[df_bivariate["cell_type"].isin(rare2_cts)].reset_index(drop = True)
        rare2_mean = df_rare2[methods_all]
        
        print (rare2_mean.shape, "==", rare2_mean.mean().values, len(rare2_mean.mean().values,))        
        grouped_df["comp_rare_def2_spatial"] = rare2_mean.mean().values
    
        num = (grouped_df["comp_rare_def1_spatial"] + grouped_df["comp_rare_def2_spatial"] + grouped_df["AUPR"])
        grouped_df["composite_rare_score"] = num/3
    else:
        grouped_df["comp_rare_def2_spatial"] = 0
        grouped_df["composite_rare_score"] = (grouped_df["comp_rare_def1_spatial"] + grouped_df["AUPR"])/2
    return grouped_df

def composite_score_cellcharter(grouped_df, df_bivariate, l_top_ct, l_bottom_ct, c_top_ct, e_top_ct, e_bottom_ct, methods_all):

    
    df_curl_top = df_bivariate[df_bivariate["cell_type"].isin(c_top_ct)][methods_all].reset_index(drop = True)
    #print (df_curl_top.mean().values, len(df_curl_top.mean().values,))
    
    grouped_df["comp_curl_top_spatial"] = df_curl_top.mean().values
    
    df_elongation_top = df_bivariate[df_bivariate["cell_type"].isin(e_top_ct)][methods_all].reset_index(drop = True)
    grouped_df["comp_elongation_top_spatial"] = df_elongation_top.mean().values
    
    df_linearity_top = df_bivariate[df_bivariate["cell_type"].isin(l_top_ct)][methods_all].reset_index(drop = True)
    grouped_df["comp_linearity_top_spatial"] = df_linearity_top.mean().values
    
    df_elongation_bottom = df_bivariate[df_bivariate["cell_type"].isin(e_bottom_ct)][methods_all].reset_index(drop = True)
    grouped_df["comp_elongation_bottom_spatial"] = df_elongation_bottom.mean().values
    
    df_linearity_bottom = df_bivariate[df_bivariate["cell_type"].isin(l_bottom_ct)][methods_all].reset_index(drop = True)
    grouped_df["comp_linearity_bottom_spatial"] = df_linearity_bottom.mean().values

    num = (grouped_df["comp_curl_top_spatial"] + grouped_df["comp_elongation_top_spatial"] + grouped_df["comp_linearity_top_spatial"] + \
    grouped_df["comp_elongation_bottom_spatial"] + grouped_df["comp_linearity_bottom_spatial"])
    
    grouped_df["composite_cellcharter_score"] = num/5

    return grouped_df

def create_label_column(grouped_df):

    new_df = pd.DataFrame()

    for col in grouped_df.columns[1:]:
        new_df[col + "_label"] = np.round(grouped_df[col].values, 2)
        
    grouped_df = pd.concat([grouped_df, new_df], axis = 1)
    
    return grouped_df

def add_aupr(grouped_df, methods_all, metrics_path):
    
    m = "AUPR"
    ### Add an index corresponding to AUPR metric from the saved files.
    data = pd.read_csv(metrics_path + m + ".csv")
    data = data[methods_all]
    data = data.dropna(how = "all")    
    data = data.mean().T.tolist() # by default mean calculation ignores the nan/none present in the array.
    print (data)
    grouped_df.loc["AUPR"] = data
    return grouped_df
    
def create_composite_score_file(metrics_path, cc_path, sheet_path, data_name, methods_all, bi_metrics, \
                                acc_metrics, spot_metrics, data_csv_path):

    all_metrics = bi_metrics + acc_metrics
    l_top_ct, l_bottom_ct, c_top_ct, e_top_ct, e_bottom_ct, rare2_cts, rare1_cts = filter_cell_types(cc_path, sheet_path, data_name)

    df_lsts = []
    for m in all_metrics:
        
        data = pd.read_csv(metrics_path + m + ".csv")
        #print (data.columns)
        data["cell_type"] = data[data.columns[0]].astype(str)
        data["Metric"] = m
        data = data[["Metric", "cell_type"] + methods_all]
        df_lsts += [data]

    df = pd.concat(df_lsts, ignore_index = True)
    
    merged_df = scaing_of_dataframe(df, methods_all, all_metrics)
    df_bivariate = merged_df[merged_df["Metric"].isin(bi_metrics)].reset_index(drop  = True) # bi_metrics
    df_accuracy = merged_df[merged_df["Metric"].isin(acc_metrics)].reset_index(drop  = True) # acc_metrics
    
    grouped_df = merged_df.groupby("Metric").mean()

    #print (grouped_df, df_bivariate, df_accuracy)
     
    grouped_df = add_aupr(grouped_df, methods_all, metrics_path)
    grouped_df = grouped_df.T  
    grouped_df["id"] = grouped_df.index.values    
    grouped_df = grouped_df[['id'] + all_metrics + ["AUPR"]]
    grouped_df = grouped_df.reset_index(drop = True)

    # Read RMSE and JS divergence.
    for m in spot_metrics:
        data = pd.read_csv(metrics_path + m + ".csv")
        grouped_df[m] = np.array(data[methods_all]).reshape(-1)
        
    grouped_df['comp_acc'] = grouped_df.loc[:, acc_metrics].mean(axis = 1)    
    grouped_df['comp_bivariate'] = grouped_df.loc[:, bi_metrics].mean(axis = 1)

    grouped_df = composite_score_rare(grouped_df, df_bivariate, rare1_cts, rare2_cts, methods_all)
    grouped_df = composite_score_cellcharter(grouped_df, df_bivariate, l_top_ct, l_bottom_ct, c_top_ct, e_top_ct, e_bottom_ct, methods_all)

    comp_score_lst = ['comp_acc', 'comp_bivariate', 'composite_rare_score', 'composite_cellcharter_score']    
    grouped_df["overall_composite_score"] = grouped_df[comp_score_lst].sum(axis = 1)/len(comp_score_lst)
    grouped_df = create_label_column(grouped_df)
    grouped_df.to_csv(data_csv_path)
    return grouped_df


def create_plotting_files(methods_all, col_group_csv_path, col_info_csv_path, row_group_csv_path, row_info_csv_path):
    
    columns_group = create_columns_group()
    columns_group.to_csv(col_group_csv_path, index = False)
    
    columns = create_columns_info(methods_all)
    columns.to_csv(col_info_csv_path, index = False)
    
    row_groups, rows_info = create_rows_group_and_info(methods_all)
    row_groups.to_csv(row_group_csv_path, index = False)
    rows_info.to_csv(row_info_csv_path, index = False)
    
    
































