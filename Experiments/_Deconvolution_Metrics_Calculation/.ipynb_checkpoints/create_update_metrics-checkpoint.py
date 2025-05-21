import pandas as pd
import numpy as np
import scanpy as sc
import scipy
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import entropy
from scipy.spatial.distance import jensenshannon
from sklearn.preprocessing import MinMaxScaler
import seaborn as sns
import os
import pickle
from matplotlib import rcParams
from scipy.spatial import distance
import sys
from metrics import *


### Evaluate all methods for a dataset...

def intersection_of_spots(gt, pred, coords):
    
    if len(gt) == len(pred):
        pass
    else:
        
        if sum(np.isin(pred.index.astype(str), gt.index.astype(str))) == 0:
            print("Error: Index not present, unable to subset")
            #sys.exit(0)
            return gt, pred
        else:
            print ("Subsetting dataset")
            spot_codes = np.intersect1d(
                np.array(pred.index.astype(str)), np.array(gt.index.astype(str))
            )
            gt = gt.loc[spot_codes]
            pred = pred.loc[spot_codes]
            coords = np.array(coords.loc[spot_codes])
            print (len(gt), len(pred), len(coords))

    return gt, pred, coords
    
def preprocess_predictions(gt, pred, method, coords):

    # Give the list of all possible column names that indicate bar codes.
    if (pred.columns[0] in ["spot_id", "Unnamed:0", "barcode"]):
        print ("bar codes as column")
        pred.index = pred[pred.columns[0]].astype(str)
    else:
        pred.index = pred.index.astype(str)
        
    pred.columns = pred.columns.astype(str)
    pred = pred.sort_index(axis = 1)
    pred.columns = pred.columns.str.replace('"', '')
    pred.columns = pred.columns.str.strip()
    
    for ch in ["-", " ", "/", ".", "+", "(", ")", ",", "&"]:        
        pred.columns = np.char.replace(np.array(pred.columns).astype(str), ch, "_")
    
    #print ("Testing column", pred.columns)

    pred = pred[~pred.isnull().any(axis=1)]
    print("Num cells dropped:" + str(len(gt) - len(pred)))
    gt, pred, coords = intersection_of_spots(gt, pred, coords)

    if method == "Cell2location":
        pred.columns = np.char.replace(
            np.array(pred.columns).astype(str), "q05cell_abundance_w_sf_", ""
        )
    if method == "SD":
        pred.drop(columns="X", inplace=True)
        
    # Making sure columns are same between ground truth and predictions.
    pred = pred[np.sort(gt.columns)]
    gt = gt[np.sort(gt.columns)]
    #print ("preprocessing in between", pred.shape, "and", gt.shape)
    
    if (method == "Polaris") or (method == "Cell2location"):
        print ("Normalizating results")
        pred = pred.div(pred.sum(axis = 1), axis = 0)

    if method == "Autogenes":
        
        pred[pred < 0] = 0
        pred = pred.div(pred.sum(axis=1), axis = 0)        
        gt = gt[np.array(~pred.isnull().any(axis = 1))]
        coords = coords[np.array(~pred.isnull().any(axis = 1))]
        pred = pred[np.array(~pred.isnull().any(axis = 1))]
    
    for col in np.array(gt.columns):
        if col in np.array(pred.columns):
            pass
        else:
            print(
                "WARNING: Cell type "
                + col
                + " not found for "
                + method
                + ", assuming all 0s"
            )
            pred[col] = np.zeros(len(pred))
       
    pred = pred[gt.columns]

    print ("preprocessing finished", pred.shape, "and", gt.shape)
    return gt, pred, coords


def evaluate_method(gt, pred, coords, method, di, dataset_name, l, eps, co, metrics):

    ### RMSE Computation; joint excel sheet not working
    """
    rmse_path = "/data/Ajita/Spatial/Datasets/Spatial_Deconvolution/_Evaluation/Metrics_Calculation/rmse.csv"
    jsd_path = "/data/Ajita/Spatial/Datasets/Spatial_Deconvolution/_Evaluation/Metrics_Calculation/jsd.csv"

    rmse = pd.read_csv(rmse_path, index_col=0)
    jsd = pd.read_csv(jsd_path, index_col=0)
    """
    gt, pred, coords = preprocess_predictions(gt, pred, method, coords)

    di["RMSE"].loc[1, method] = get_rmse(pred, gt)
    di["JSD"].loc[1, method] = get_jsd(pred, gt)
    """
    rmse.loc[dataset_name, method] = di["RMSE"].loc[1, method]
    jsd.loc[dataset_name, method] = di["JSD"].loc[1, method]
    rmse.to_csv(rmse_path)
    jsd.to_csv(jsd_path)
    """
    # Other metrics computation

    # Check to see if only RMSE and JS needs to be updated.By default, RMSE and JS will always be updated.
    if (len(metrics) <= 2) and (("RMSE" in metrics) or ("JSD" in metrics)):
        print ("only RMSE and JS updated")
        return di
        
    for col in np.array(gt.columns):

        for m in metrics: 
            
            if (m == "pearson"):
                val = get_pearson(gt[col], pred[col], eps = eps)
                di["pearson"].loc[col, method] = val
                #if (method == "STRIDE"):
                #    print ("pearson", val)
            
            elif (m == "cosine_sim"):
                val = get_cosine_sim(gt[col], pred[col], eps)
                di["cosine_sim"].loc[col, method] = val
                #print ("cosine similarity", val)

            elif (m == "morans_r"):
                val = get_morans_R(gt[col], pred[col], coords, l=l, co=co, eps = eps)
                di["morans_r"].loc[col, method] = val
                #if (method == "STRIDE"):
                #print ("Moran's R", val)
                # di['spearman'].loc[col, method] = get_spearman(gt[col],pred[col],coords,l = l, co = co)

            elif(m == "spatial_pearson"):
                val = get_spatial_pearson(gt[col], pred[col], coords, l=l, co=co, eps = eps)
                di["spatial_pearson"].loc[col, method] = val
                #print ("spatial pearson", val)
            
            elif (m == "ssim"):
                val = compute_ssim(gt[col], pred[col], eps = eps)
                di["ssim"].loc[col, method] = val
                #print ("ssim", val)

            elif (m == "lee_stat"):
                val = compute_Lee_stats(gt[col], pred[col], coords, l=l, co=co, eps = eps)
                di["lee_stat"].loc[col, method] = val
                #print ("lee stats", val)

            elif (m == "geary_c"):
                val = compute_geary(gt[col], pred[col], coords, l=l, co=co, eps = eps)
                di["geary_c"].loc[col, method] = val
                #print ("geary_c", val)
                
            elif (m == "AUPR"):
                compute_AUPR(gt[col], pred[col], coords, l=l, co=co, eps = eps)
    
    print (method + "done")
    return di

def preprocess_groundtruth(ground_truth_path, col, dataset_name):
    
    adata = sc.read_h5ad(ground_truth_path)

    if adata.var_names[0][0:2] == "q0":

        tmp = []
        for k in range(len(adata.var_names)):
            tmp.append(adata.var_names[k][23:])
        adata.var_names = tmp

    for ch in ["-", " ", "/", ".", "+", "(", ")", ",", "&"]:
        adata.var_names = adata.var_names.str.strip()
        adata.var_names =  adata.var_names.str.replace('"', '')
        adata.var_names = np.char.replace(
            np.array(adata.var_names).astype(str), ch, "_"
        )

    if (isinstance(adata.obsm[col], pd.DataFrame)):
        gt = adata.obsm[col]
        gt.columns = adata.var_names
        gt.index = adata.obs_names
        
    else:
        gt = pd.DataFrame(adata.obsm[col], columns = adata.var_names, index = adata.obs_names)
        
    gt.columns = gt.columns.astype(str)

    # coords...
    if "array_row" in adata.obs.columns:
        coords = adata.obs[["array_row", "array_col"]]  
    elif "new_x" in adata.obs.columns:
        coords = adata.obs[["new_x", "new_y"]]
    elif "xcoord" in adata.obs.columns:
        coords = adata.obs[["xcoord", "ycoord"]]
    elif "x_new" in adata.obs.columns:
        coords = adata.obs[["x_new", "y_new"]]
    elif "Centroid_X" in adata.obs.columns:
        coords = adata.obs[["Centroid_X", "Centroid_Y"]]
    elif "center_x" in adata.obs.columns:
        coords = adata.obs[["center_x", "center_y"]]
    elif "x" in adata.obs.columns:
        coords = adata.obs[["x", "y"]]

    dataset_lst = ["Merfish_brain", "Merfish_ileum", "Merfish_Lung_Cancer", "Merfish_Breast_Cancer"]
    if (("Breast_Atlas" in dataset_name) or ("Kidney_Atlas" in dataset_name)) or (dataset_name in dataset_lst):
        
        print ("Using pixel coordinates", dataset_name)
        coords = pd.DataFrame(
            [adata.obsm["spatial"][:, 1], adata.obsm["spatial"][:, 0]]
        ).T
        coords.columns = ["new_x", "new_y"]
        coords.index = adata.obs.index
        
    """
    # I want code to give error when spot coordinates are not present.
    else:
        coords = pd.DataFrame(
            [adata.obsm["spatial"][:, 1], adata.obsm["spatial"][:, 0]]
        ).T
        coords.columns = ["new_x", "new_y"]
        coords.index = adata.obs.index
    """        
    return gt, coords, adata


def create_new_evaluation(path_to_outputs, ground_truth_path, rmse_path, jsd_path, dataset_name, col, methods_all, celltype_metrics, 
                          global_metrics, l = 1.2, eps = 1e-8, co = 0):
    
    os.makedirs(path_to_outputs + 'Plots/', exist_ok = True)
    os.makedirs(path_to_outputs + 'Metrics/', exist_ok = True)

    """
    rmse = pd.read_csv(rmse_path, index_col = 0)
    jsd = pd.read_csv(jsd_path, index_col = 0)
    rmse.loc[dataset_name] = np.zeros(len(methods_all))
    jsd.loc[dataset_name] = np.zeros(len(methods_all))
    rmse.to_csv(rmse_path)
    jsd.to_csv(jsd_path)
    """
    gt, coords, adata = preprocess_groundtruth(ground_truth_path, col, dataset_name)

    random_matrix = pd.DataFrame(np.random.rand(gt.shape[0],gt.shape[1]), columns = gt.columns, index = gt.index)
    random_matrix = random_matrix.div(random_matrix.sum(axis=1), axis=0)
    random_matrix.to_csv(path_to_outputs + "output_zRandom.csv")
    
    di = {}
    for metric in (celltype_metrics):
        di[metric] = pd.DataFrame(columns = methods_all, index = adata.var_names)
    for metric in (global_metrics):
        di[metric] = pd.DataFrame(columns = methods_all, index = [1])
        
    for method in methods_all:
        print (method)
        try:
            if method == "STRIDE":
                pred = pd.read_table(
                    path_to_outputs + "output_" + method + ".csv", index_col=0, sep="\t"
                )
            elif method == "Polaris":
                pred = pd.read_table(
                    path_to_outputs + "output_" + method + ".tsv", index_col=0
                )
            else:
                pred = pd.read_csv(
                    path_to_outputs + "output_" + method + ".csv", index_col=0
                )
        except:
            print(method + " output not found")
            continue
        # If new method is getting evaluated, all metrics should be updated.
        di = evaluate_method(gt,pred,coords,method,di, dataset_name, l, eps, co, celltype_metrics + global_metrics)
        
    pickle.dump(di,open(path_to_outputs + "Metrics/eval.pkl","wb"))
    
    for metric in (celltype_metrics + global_metrics):
        di[metric].to_csv(path_to_outputs + "Metrics/" + metric + '.csv')


def update_method_evaluation(path_to_outputs, ground_truth_path, rmse_path, jsd_path, method, dataset_name, col, celltype_metrics, global_metrics, l, eps, co):

    gt, coords, adata = preprocess_groundtruth(ground_truth_path, col, dataset_name)
    di = pickle.load(open(path_to_outputs + "Metrics/eval.pkl", "rb"))

    print (method)
    try:
        if method == "STRIDE":
            pred = pd.read_table(
                path_to_outputs + "output_" + method + ".csv", index_col=0, sep="\t"
            )
        elif method == "Polaris":
            pred = pd.read_table(
                path_to_outputs + "output_" + method + ".tsv", index_col=0
            )
        else:
            pred = pd.read_csv(
                path_to_outputs + "output_" + method + ".csv", index_col=0
            )
    except:
        print(method + " output not found")
        return
            
    di = evaluate_method(gt, pred, coords, method, di, dataset_name, l, eps, co, celltype_metrics)
        
    pickle.dump(di, open(path_to_outputs + "Metrics/eval.pkl", "wb"))
    
    for metric in (celltype_metrics + global_metrics):
        if (metric != "AUPR"):
            di[metric].to_csv(path_to_outputs + "Metrics/" + metric + ".csv")
            
def update_metric_evaluation(path_to_outputs, ground_truth_path, rmse_path, jsd_path, dataset_name, col, methods_all, metrics, global_metrics, l, eps, co):

    gt, coords, adata = preprocess_groundtruth(ground_truth_path, col, dataset_name)
    di = pickle.load(open(path_to_outputs + "Metrics/eval.pkl", "rb"))

    for method in methods_all:

        print (method)
        try:
            if method == "STRIDE":
                pred = pd.read_table(
                    path_to_outputs + "output_" + method + ".csv", index_col=0, sep="\t"
                )
            elif method == "Polaris":
                pred = pd.read_table(
                    path_to_outputs + "output_" + method + ".tsv", index_col=0
                )
            else:
                pred = pd.read_csv(
                    path_to_outputs + "output_" + method + ".csv", index_col=0
                )
        except:
            print(method + " output not found")
            continue
            
        di = evaluate_method(gt, pred, coords, method, di, dataset_name, l, eps, co, metrics)
        
    pickle.dump(di, open(path_to_outputs + "Metrics/eval.pkl", "wb"))
    
    for metric in metrics:
        if (metric != "AUPR"):
            di[metric].to_csv(path_to_outputs + "Metrics/" + metric + ".csv")
        

def plot_histogram(ground_truth_path, dataset_name, col, l, eps, co, rare1_cts, rare2_cts):

    gt, coords, adata = preprocess_groundtruth(ground_truth_path, col, dataset_name)

    print (gt.columns)
    
def make_plots(path_to_outputs, dataset_name, celltype_metrics):
    
    rcParams['figure.figsize'] = 20,8
    di = pickle.load(open(path_to_outputs + 'Metrics/eval.pkl',"rb"))
    
    for metric in celltype_metrics:
        
        di[metric] = di[metric].dropna(axis=1, how='all')
        ymin = np.min(np.array(di[metric]))
        ymax = np.max(np.array(di[metric]))
        ymax = 1
        ymin = 0
        ax = sns.boxplot(data = di[metric])
        ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
        plt.ylim(ymin-0.1*ymax,ymax + 0.1*ymax)
        plt.xlabel('Method')
        plt.ylabel(metric)
        plt.title(dataset_name + "_" + metric)
        plt.savefig(path_to_outputs + "Plots/" +metric+".png")
        plt.show()

        