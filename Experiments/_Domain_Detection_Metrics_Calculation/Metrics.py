# ----------------------- Imports -----------------------


import os
import re
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import anndata as an
import seaborn as sns
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib import font_manager as fm
from scipy.spatial import distance_matrix
from scipy.spatial.distance import squareform, pdist
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
    homogeneity_score,
    completeness_score
)
from sklearn.preprocessing import StandardScaler

# ----------------------- Utility -----------------------
def fx_1NN(i,location_in):
        location_in = np.array(location_in)
        dist_array = distance_matrix(location_in[i,:][None,:],location_in)[0,:]
        dist_array[i] = np.inf
        return np.min(dist_array)
    

def fx_kNN(i,location_in,k,cluster_in):

    location_in = np.array(location_in)
    cluster_in = np.array(cluster_in)


    dist_array = distance_matrix(location_in[i,:][None,:],location_in)[0,:]
    dist_array[i] = np.inf
    ind = np.argsort(dist_array)[:k]
    cluster_use = np.array(cluster_in)
    if np.sum(cluster_use[ind]!=cluster_in[i])>(k/2):
        return 0
    else:
        return 1

def _compute_CHAOS(clusterlabel, location):

    clusterlabel = np.array(clusterlabel)
    location = np.array(location)
    matched_location = StandardScaler().fit_transform(location)

    clusterlabel_unique = np.unique(clusterlabel)
    dist_val = np.zeros(len(clusterlabel_unique))
    count = 0
    for k in clusterlabel_unique:
        location_cluster = matched_location[clusterlabel==k,:]
        if len(location_cluster)<=2:
            continue
        n_location_cluster = len(location_cluster)
        results = [fx_1NN(i,location_cluster) for i in range(n_location_cluster)]
        dist_val[count] = np.sum(results)
        count = count + 1

    chaos = np.sum(dist_val)/len(clusterlabel)
    return np.exp(-0.5*chaos)

def _compute_PAS(clusterlabel,location):

    clusterlabel = np.array(clusterlabel)
    location = np.array(location)
    matched_location = location
    results = [fx_kNN(i,matched_location,k=10,cluster_in=clusterlabel) for i in range(matched_location.shape[0])]
    return np.sum(results)/len(clusterlabel)

def compute_ASW(pred,spatial_coords):
    distance_matrix = squareform(pdist(spatial_coords))
    sil =  silhouette_score(X=distance_matrix, labels=pred, metric='precomputed')
    return (sil + 1)/2
    
def LISI(coords, meta, label, perplexity=30, nn_eps=0):
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    from rpy2.robjects.packages import importr
    importr("lisi")
    if not isinstance(coords, pd.DataFrame):
        coords = pd.DataFrame(coords)
    if not isinstance(meta, pd.DataFrame):
        meta = pd.DataFrame(meta)
    meta = meta.loc[:, [label]]
    meta[label] = meta[label].astype(str)

    coords = robjects.conversion.py2rpy(coords)
    meta = robjects.conversion.py2rpy(meta)
    as_matrix = robjects.r["as.matrix"]
    lisi = robjects.r["compute_lisi"](as_matrix(coords), meta, label, perplexity, nn_eps)
    if isinstance(lisi, pd.DataFrame):
        lisi = lisi.values
    elif isinstance(lisi, np.recarray):
        lisi = [item[0] for item in lisi]

    return lisi

# ----------------------- Import functions -----------------------
def import_dataset(dataset_name, paths, mode = "evaluation"):
    """
    Load dataset given its name and a dictionary of file paths.
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset (used for branching logic).
    paths : dict
        Dictionary of file paths needed for that dataset.
        Must contain at least 'st_path'.
        Can contain 'gnd_path' if ground truth is from a CSV/TSV file.
    
    Returns
    -------
    gnd : pd.Series or pd.DataFrame
        Ground truth clusters.
    locs : pd.DataFrame
        Spatial coordinates.
    """
    # Load main AnnData object
    ST = an.read_h5ad(paths["st_path"])

    if dataset_name.startswith("DLPFC"):
        gnd = pd.read_csv(paths["gnd_path"], sep="\t")['layer_guess_reordered']
        gnd.index.name = "Index"
        locs = pd.DataFrame({
            "array_row": ST.obs.array_row,
            "array_col": ST.obs.array_col
        })

    elif dataset_name.startswith("embryo"):
        gnd = pd.DataFrame(ST.obs['annotation']).rename(columns={"annotation": "Cluster"})
        gnd.index.name = "Index"
        locs = pd.DataFrame({
            "array_row": ST.obsm["spatial"][:, 0],
            "array_col": ST.obsm["spatial"][:, 1]
        }, index=ST.obs_names)

    elif dataset_name == "mouse_breast_cancer":
        gnd = pd.DataFrame(ST.obs['ground_truth']).rename(columns={"annotation": "Cluster"})
        gnd.index.name = "Index"
        locs = pd.DataFrame({
            "array_row": ST.obs.array_row,
            "array_col": ST.obs.array_col
        })

    elif dataset_name.startswith("MERFISH_brain"):
        gnd = pd.DataFrame(ST.obs['ground_truth']).rename(columns={"annotation": "Cluster"})
        gnd.index.name = "Index"
        locs = pd.DataFrame({
            "array_row": ST.obsm["spatial"][:, 0],
            "array_col": ST.obsm["spatial"][:, 1]
        }, index=ST.obs_names)

    elif dataset_name == "simulated_breast_atlas":
        gnd = pd.read_csv(paths["gnd_path"]).set_index(ST.obs_names)['Ground Truth']
        gnd.index.name = "Index"
        locs = pd.DataFrame({
            "array_row": ST.obsm["spatial"][:, 0],
            "array_col": ST.obsm["spatial"][:, 1]
        }, index=ST.obs_names)

    elif dataset_name == "osmFISH":
        gnd = pd.DataFrame(ST.obs['ground_truth']).rename(columns={"annotation": "Cluster"})
        gnd.index.name = "Index"
        locs = pd.DataFrame({
            "array_row": ST.obsm["spatial"][:, 0],
            "array_col": ST.obsm["spatial"][:, 1]
        }, index=ST.obs_names)

    elif dataset_name.startswith("simulated_kidney_cancer"):
        gnd = pd.read_csv(paths["gnd_path"]).set_index(ST.obs_names)['Ground Truth']
        gnd.index.name = "Index"
        locs = pd.DataFrame({
            "array_row": ST.obs.new_x,
            "array_col": ST.obs.new_y
        })

    elif dataset_name.startswith("simulated_breast_cancer"):
        ST.obs_names_make_unique()
        gnd = pd.read_csv(paths["gnd_path"]).set_index(ST.obs_names)['Ground Truth']
        gnd.index.name = "Index"
        locs = pd.DataFrame({
            "array_row": ST.obs.new_x,
            "array_col": ST.obs.new_y
        })

    elif dataset_name.startswith("simulated_liver_cancer"):
        gnd = pd.read_csv(paths["gnd_path"]).set_index(ST.obs_names)['Ground Truth']
        gnd.index.name = "Index"
        locs = pd.DataFrame({
            "array_row": ST.obs.new_x,
            "array_col": ST.obs.new_y
        })

    elif dataset_name.startswith("simulated_intestine"):
        gnd = pd.read_csv(paths["gnd_path"]).set_index(ST.obs_names)['Ground Truth']
        gnd.index.name = "Index"
        locs = pd.DataFrame({
            "array_row": ST.obs.new_x,
            "array_col": ST.obs.new_y
        })

    elif dataset_name == "simulated_chicken_heart":
        gnd = pd.read_csv(paths["gnd_path"]).set_index(ST.obs_names)['Ground Truth']
        gnd.index.name = "Index"
        locs = pd.DataFrame({
            "array_row": ST.obs.array_row,
            "array_col": ST.obs.array_col
        })

    elif dataset_name == "simulated_prostate_cancer":
        gnd = pd.read_csv(paths["gnd_path"]).set_index(ST.obs_names)['Ground Truth']
        gnd.index.name = "Index"
        locs = pd.DataFrame({
            "array_row": ST.obs.x_new,
            "array_col": ST.obs.y_new
        })

    elif dataset_name == "simulated_cerebellum":
        gnd = pd.read_csv(paths["gnd_path"]).set_index(ST.obs_names.astype(int))['Ground Truth']
        gnd.index.name = "Index"
        locs = pd.DataFrame({
            "array_row": ST.obs.xcoord,
            "array_col": ST.obs.ycoord
        }).set_index(ST.obs_names.astype(int))

    else:
        raise ValueError(f"Incorrect input: {dataset_name}")
    
    if mode == "plot":
        locs = pd.DataFrame({
            "array_row": ST.obsm["spatial"][:, 0],
            "array_col": ST.obsm["spatial"][:, 1]
        }).set_index(ST.obs_names)
        if dataset_name == "simulated_cerebellum":
            locs.index = locs.index.astype(int)
    locs.index.name = "Index"
    
    return gnd, locs


def load_prediction(dataset_name, method_name, path):
    """
    Load output file given its name and file path.
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset (used for branching logic).
    method_name : str
        Name of the method (used for branching logic).
    path : str
        Path to output file
    Returns
    -------
    pred : predicted clusters
    """

    if os.path.exists(path):

        if dataset_name.startswith("simulated_breast_cancer") and method_name in ["BASS","PRECAST", "BayesSpace" , "DR_SC"]:
            pred = pd.read_csv(path)
            pred = pred.set_index(pred.columns[0])
            if method_name == "PRECAST":
                pred = pred['cluster']
            pred.index = pred.index.str.replace(r'_\d+$', '', regex=True)

        elif(method_name == "banksy"):
            pred = pd.read_csv(path)
            pred = pred.set_index(pred.columns[0])
            columns_to_grab = [col for col in pred.columns if col.startswith("clust_M1_lam0.")]
            pred = pred[columns_to_grab]

        elif method_name == "giotto":     
            pred = pd.read_csv(path)
            if "leiden_clus" in pred.columns:
                columns_to_grab = ["leiden_clus"]
                pred = pred.set_index(pred.columns[1])
            else:
                columns_to_grab = ["cluster"]
                pred = pred.set_index(pred.columns[0])
            pred = pred[columns_to_grab]
            
        elif method_name == "PRECAST" or method_name == "BayesCafe":
            pred = pd.read_csv(path)
            pred = pred.set_index(pred.columns[0])
            columns_to_grab = ["cluster"]
            pred = pred[columns_to_grab]

        elif method_name == "SpaceFlow":
            pred = pd.read_csv(path)
            pred = pred.set_index(pred.columns[1])
            columns_to_grab = ["Predicted_cell_label"]
            pred = pred[columns_to_grab]

        else:
            pred = pd.read_csv(path)
            pred = pred.set_index(pred.columns[0])

    else:
        print(f"path doesn't exist : {path}")
        print(f"{method_name} output missing for {dataset_name}")
        pred = None

    return pred


# ----------------------- Metric Computationn -----------------------


def compute_metrics(dataset_names, dataset_paths, method_names, pred_paths, error_log_file, output_dir, metric_names="all"):
    """
    Compute metric values given dataset and method names and paths
    
    Parameters
    ----------
    dataset_names : list
        List of datasets to compute on
        Options for dataset names present at bottom of this file
    dataset_paths: dictionary of dictionaries
        Dictionaries mapping dataset names to dictionaries with st_path and gnd_path needed for importing dataset
    method_names: list
        List of methods to compute on
        Options for method names present at bottom of this file
    pred_paths : dictionary of dictionaries
        Dictionary mapping dataset names and method names to path of output file
    metric_names : list or str
        List of metrics to compute or "all"
    error_log_file : str
        Path to error log file
    output_dir : str
        Path to save final metric compute files to
    Returns
    -------
    pred : predicted clusters
    """
    # Available metrics
    available_metrics = ["ARI", "NMI", "CHAOS", "PAS", "ASW", "HOM", "COM"]
    if metric_names == "all":
        metric_names = available_metrics
    
    # Initialize pivoted DataFrame for each metric
    method_names_dict = {
    "SCANIT": "SCANIT",
    "CCST": "CCST",
    "DeepST": "DeepST",
    "GraphST": "GraphST",
    "PROST": "PROST",
    "SpaSRL": "SpaSRL",
    "STAGATE": "STAGATE",
    "SpatialPCA": "SpatialPCA",
    "banksy": "Banksy",
    "giotto": "Giotto",
    "DR_SC": "DR.SC",
    "ISC_MEB": "ISC.MEB",
    "BayesSpace": "BayesSpace",
    "PRECAST": "PRECAST",
    "BayesCafe": "BayesCafe",
    "BASS": "BASS",
    "SpaceFlow" : "SpaceFlow",
    "IRIS" : "IRIS"
    }

    dataset_names_dict = {
    "DLPFC151507": "DLPFC 151507",
    "DLPFC151508": "DLPFC 151508",
    "DLPFC151509": "DLPFC 151509",
    "DLPFC151510": "DLPFC 151510",
    "DLPFC151669": "DLPFC 151669",
    "DLPFC151670": "DLPFC 151670",
    "DLPFC151671": "DLPFC 151671",
    "DLPFC151672": "DLPFC 151672",
    "DLPFC151673": "DLPFC 151673",
    "DLPFC151674": "DLPFC 151674",
    "DLPFC151675": "DLPFC 151675",
    "DLPFC151676": "DLPFC 151676",
    "embryo9.5": "Embryo 9.5",
    "embryo14.5": "Embryo 14.5",
    "mouse_breast_cancer": "Mouse Breast Cancer",
    "MERFISH_brain0.04": "MERFISH Brain 0.04",
    "MERFISH_brain0.09": "MERFISH Brain 0.09",
    "MERFISH_brain0.14": "MERFISH Brain 0.14",
    "MERFISH_brain0.19": "MERFISH Brain 0.19",
    "MERFISH_brain0.24": "MERFISH Brain 0.24",
    "osmFISH": "osmFISH",
    "simulated_kidney_cancer410": "Kidney Cancer 410",
    "simulated_kidney_cancer411": "Kidney Cancer 411",
    "simulated_kidney_cancer506": "Kidney Cancer 506",
    "simulated_breast_cancerER+_CID4290": "Breast Cancer ER+ CID 4290",
    "simulated_breast_cancerTNBC_CID44971": "Breast Cancer TNBC CID 44971",
    "simulated_liver_cancerHCC-1L": "Liver Cancer HCC-1L",
    "simulated_liver_cancerHCC-2L": "Liver Cancer HCC-2L",
    "simulated_liver_cancerHCC-3L": "Liver Cancer HCC-3L",
    "simulated_liver_cancerHCC-4L": "Liver Cancer HCC-4L",
    "simulated_breast_atlas": "Breast Atlas",
    "simulated_intestineA1": "Intestine A1",
    "simulated_intestineA2": "Intestine A2",
    "simulated_chicken_heart": "Chicken Heart",
    "simulated_prostate_cancer" : "Prostate Cancer",
    "simulated_cerebellum" : "Cerebellum"
    }



    metric_results = {
        metric: pd.DataFrame(
            index=list(method_names_dict.values()), columns=list(dataset_names_dict.values()), dtype=float
        ) for metric in metric_names
    }

    os.makedirs(output_dir, exist_ok=True)
    print(f"saving to {output_dir}")

    # Process each dataset
    for dataset_name in dataset_names:
        try:
            # Load dataset (ground truth and spatial coordinates)
            gnd, locs = import_dataset(dataset_name , dataset_paths[dataset_name])
            gnd = gnd.dropna()
        except Exception as e:
            with open(error_log_file, "a") as log:
                log.write(f"Dataset loading error: {dataset_name} - {e}\n")
            continue

        for method_name in method_names:
            try:
                # Load predictions
                print(f"Running on {dataset_name}_{method_name}")
                pred = load_prediction(dataset_name, method_name, pred_paths[dataset_name][method_name])

                # If predictions are missing, skip to the next method
                if pred is None:
                    continue
                pred = pred[~pred.index.duplicated(keep='first')]

                # Ensure first columns are indices
                intersect_idx = gnd.index.intersection(pred.index)

                # Filter gnd and pred based on the intersection
                gnd_filtered = gnd.loc[intersect_idx]
                pred_filtered = pred.loc[intersect_idx]

                # Replace pred and gnd for further computation
                gnd_values = gnd_filtered.values.flatten()
                pred_values = pred_filtered.iloc[:, 0].values

                # Placeholder for spatial coordinates
                spatial_coords = locs.loc[intersect_idx].values

                # Compute metrics and populate the pivoted DataFrame
                for metric in metric_names:
                    try:
                        if metric == "ARI":
                            value = adjusted_rand_score(gnd_values, pred_values)
                        elif metric == "NMI":
                            value = normalized_mutual_info_score(gnd_values, pred_values)
                        elif metric == "CHAOS":
                            value = _compute_CHAOS(pred_values, spatial_coords)
                        elif metric == "PAS":
                            value = _compute_PAS(pred_values, spatial_coords)
                        elif metric == "ASW":
                            value = compute_ASW(pred_values, spatial_coords)
                        elif metric == "HOM":
                            value = homogeneity_score(gnd_values, pred_values)
                        elif metric == "COM":
                            value = completeness_score(gnd_values, pred_values)
                        else:
                            continue

                        # Populate the DataFrame
                        metric_results[metric].loc[method_names_dict[method_name], dataset_names_dict[dataset_name]] = value

                    except Exception as e:
                        with open(error_log_file, "a") as log:
                            log.write(f"Metric error: Dataset={dataset_name}, Method={method_name}, Metric={metric} - {e}\n")
            except Exception as e:
                with open(error_log_file, "a") as log:
                    log.write(f"Prediction loading error: Dataset={dataset_name}, Method={method_name} - {e}\n")

    # Save each metric result as a CSV
    for metric, df in metric_results.items():
        try:
            output_file = os.path.join(output_dir, f"{metric}_results.csv")
            df.to_csv(output_file)
            print(f"Saved {metric} results to {output_file}")
        except Exception as e:
            with open(error_log_file, "a") as log:
                log.write(f"Metric saving error: Metric={metric} - {e}\n")

def update_metrics(dataset_names, dataset_paths, method_names, pred_paths, error_log_file , output_dir,metric_names="all"):

    """
    Modify existing output files with new metric values given dataset and method names and paths
    
    Parameters
    ----------
    dataset_names : list
        List of datasets to compute on
        Options for dataset names present at bottom of this file
    dataset_paths: dictionary of dictionaries
        Dictionaries mapping dataset names to dictionaries with st_path and gnd_path needed for importing dataset
    method_names: list
        List of methods to compute on
        Options for method names present at bottom of this file
    pred_paths : dictionary of dictionaries
        Dictionary mapping dataset names and method names to path of output file
    metric_names : list or str
        List of metrics to compute or "all"
    error_log_file : str
        Path to error log file
    output_dir : str
        Path to save final metric compute files to
    Returns
    -------
    pred : predicted clusters
    """
    # Available metrics
    available_metrics = ["ARI", "NMI", "CHAOS", "PAS", "ASW", "HOM", "COM"]
    if metric_names == "all":
        metric_names = available_metrics
    method_names_dict = {
    "SCANIT": "SCANIT",
    "CCST": "CCST",
    "DeepST": "DeepST",
    "GraphST": "GraphST",
    "PROST": "PROST",
    "SpaSRL": "SpaSRL",
    "STAGATE": "STAGATE",
    "SpatialPCA": "SpatialPCA",
    "banksy": "Banksy",
    "giotto": "Giotto",
    "DR_SC": "DR.SC",
    "ISC_MEB": "ISC.MEB",
    "BayesSpace": "BayesSpace",
    "PRECAST": "PRECAST",
    "BayesCafe": "BayesCafe",
    "BASS": "BASS",
    "SpaceFlow":"SpaceFlow",
    "IRIS":"IRIS"
    }

    dataset_names_dict = {
    "DLPFC151507": "DLPFC 151507",
    "DLPFC151508": "DLPFC 151508",
    "DLPFC151509": "DLPFC 151509",
    "DLPFC151510": "DLPFC 151510",
    "DLPFC151669": "DLPFC 151669",
    "DLPFC151670": "DLPFC 151670",
    "DLPFC151671": "DLPFC 151671",
    "DLPFC151672": "DLPFC 151672",
    "DLPFC151673": "DLPFC 151673",
    "DLPFC151674": "DLPFC 151674",
    "DLPFC151675": "DLPFC 151675",
    "DLPFC151676": "DLPFC 151676",
    "embryo9.5": "Embryo 9.5",
    "embryo14.5": "Embryo 14.5",
    "mouse_breast_cancer": "Mouse Breast Cancer",
    "MERFISH_brain0.04": "MERFISH Brain 0.04",
    "MERFISH_brain0.09": "MERFISH Brain 0.09",
    "MERFISH_brain0.14": "MERFISH Brain 0.14",
    "MERFISH_brain0.19": "MERFISH Brain 0.19",
    "MERFISH_brain0.24": "MERFISH Brain 0.24",
    "osmFISH": "osmFISH",
    "simulated_kidney_cancer410": "Kidney Cancer 410",
    "simulated_kidney_cancer411": "Kidney Cancer 411",
    "simulated_kidney_cancer506": "Kidney Cancer 506",
    "simulated_breast_cancerER+_CID4290": "Breast Cancer ER+ CID 4290",
    "simulated_breast_cancerTNBC_CID44971": "Breast Cancer TNBC CID 44971",
    "simulated_liver_cancerHCC-1L": "Liver Cancer HCC-1L",
    "simulated_liver_cancerHCC-2L": "Liver Cancer HCC-2L",
    "simulated_liver_cancerHCC-3L": "Liver Cancer HCC-3L",
    "simulated_liver_cancerHCC-4L": "Liver Cancer HCC-4L",
    "simulated_breast_atlas": "Breast Atlas",
    "simulated_intestineA1": "Intestine A1",
    "simulated_intestineA2": "Intestine A2",
    "simulated_chicken_heart": "Chicken Heart",
    "simulated_prostate_cancer" : "Prostate Cancer",
    "simulated_cerebellum" : "Cerebellum"
    }
    # Initialize pivoted DataFrame for each metric
    metric_results = {}
    for metric in metric_names:
        path = f"{output_dir}/{metric}_results.csv"
        try:
            df = pd.read_csv(path,skip_blank_lines=True)
            df.set_index(df.columns[0], inplace=True, drop=True)  # Set the first column as index
            df.index.name = None  # Optional: remove index name
            metric_results[metric] = df
        except FileNotFoundError:
            # Create a new DataFrame if not already present
            metric_results[metric] = pd.DataFrame(index=method_names_dict.values(), columns=dataset_names_dict.values())
            print(f"Created new DataFrame for metric {metric}")
        

    os.makedirs(output_dir, exist_ok=True)

    # Process each dataset
    for dataset_name in dataset_names:
        try:
            # Load dataset (ground truth and spatial coordinates)
            gnd, locs = import_dataset(dataset_name , dataset_paths[dataset_name])
            gnd = gnd.dropna()

        except Exception as e:
            with open(error_log_file, "a") as log:
                log.write(f"Dataset loading error: {dataset_name} - {e}\n")
            continue

        for method_name in method_names:
            try:

                method_display_name = method_names_dict[method_name]
                # Load predictions
                print(f"Running on {dataset_name}_{method_name}")
                pred = load_prediction(dataset_name, method_name , pred_paths[dataset_name][method_name])
                

                # If predictions are missing, skip to the next method
                if pred is None:
                    continue
                pred = pred[~pred.index.duplicated(keep='first')]

                # Ensure first columns are indices
                intersect_idx = gnd.index.intersection(pred.index)

                # Filter gnd and pred based on the intersection
                gnd_filtered = gnd.loc[intersect_idx]
                pred_filtered = pred.loc[intersect_idx]

                # Replace pred and gnd for further computation
                gnd_values = gnd_filtered.values.flatten()
                
                pred_values = pred_filtered.iloc[:,0].values
                # Placeholder for spatial coordinates
                spatial_coords = locs.loc[intersect_idx].values

                # Compute metrics and populate the pivoted DataFrame
                for metric in metric_names:
                    try:
                        if metric == "ARI":
                            value = adjusted_rand_score(gnd_values, pred_values)
                        elif metric == "NMI":
                            value = normalized_mutual_info_score(gnd_values, pred_values)
                        elif metric == "CHAOS":
                            value = _compute_CHAOS(pred_values, spatial_coords)
                        elif metric == "PAS":
                            value = _compute_PAS(pred_values, spatial_coords)
                        elif metric == "ASW":
                            value = compute_ASW(pred_values, spatial_coords)
                        elif metric == "HOM":
                            value = homogeneity_score(gnd_values, pred_values)
                        elif metric == "COM":
                            value = completeness_score(gnd_values, pred_values)
                        else:
                            continue

                        # Populate the DataFrame
                        if method_display_name not in metric_results[metric].index:
                            metric_results[metric].loc[method_display_name] = np.nan
                        metric_results[metric].loc[method_names_dict[method_name], dataset_names_dict[dataset_name]] = value
                        print(f"set value of {metric} to {value}")

                    except Exception as e:
                        with open(error_log_file, "a") as log:
                            log.write(f"Metric error: Dataset={dataset_name}, Method={method_name}, Metric={metric} - {e}\n")

            except Exception as e:
                with open(error_log_file, "a") as log:
                    log.write(f"Prediction loading error: Dataset={dataset_name}, Method={method_name} - {e}\n")

    # Save each metric result as a CSV
    
    for metric, df in metric_results.items():
        try:
            output_file = os.path.join(output_dir, f"{metric}_results.csv")
            df = df.iloc[:len(method_names_dict)]  # Dynamic cutoff (good)
            df.to_csv(output_file)
            print(f"Saved {metric} results to {output_file}")
        except Exception as e:
            with open(error_log_file, "a") as log:
                log.write(f"Metric saving error: Metric={metric} - {e}\n")

# ----------------------- Visualization Functions -----------------------


def plot_stitched_heatmaps(metric_files, dataset_names, dataset_type, output_file):

    """
    Produce stitched heatmaps for all metrics
    
    Parameters
    ----------
    metric_files : dictionary
        Mapping of metric names to output files 
    dataset_names: list
        Names of datasets to include in heatmap
    dataset_type: str
        Represents name of selected group of datasets (for example : simulated/real/DLPFC)
        Used only for naming
    output_file : str
        Path to save output image to
    """
    
    num_metrics = len(metric_files) + 1
    individual_plot_height = 8  # Height of each individual heatmap
    total_height = individual_plot_height * num_metrics  # Total height for all metrics
    fig, axes = plt.subplots(4, 2, figsize=(10, total_height))  # Adjust figure size dynamically
    axes = axes.flatten()
    
    for ax, (metric, file_path) in zip(axes, metric_files.items()):
        try:
            data = pd.read_csv(file_path, index_col=0)
            data_subset = data[dataset_names]

            sns.heatmap(
                data_subset,
                cmap="coolwarm",
                cbar_kws={'label': metric},
                annot=True,  # Add annotations
                fmt=".2f",   # Format for annotations
                annot_kws={"size": 6 ,"weight": "bold" },  # Small font size for annotations
                ax=ax
            )
            ax.set_title(f"{dataset_type} Dataset Heatmap - {metric}" , fontweight = "bold")
            ax.set_xlabel("Datasets" , fontweight = "bold")
            ax.set_ylabel("Methods" , fontweight = "bold")

            labels = ax.get_xticklabels()
            bold_font = fm.FontProperties(weight='bold')
            ax.set_xticklabels(labels, rotation=45, ha="right", fontproperties=bold_font)

            y_labels = ax.get_yticklabels()
            ax.set_yticklabels(y_labels, fontproperties=bold_font)

        except Exception as e:
            ax.set_visible(False)  # Hide axes if there's an error
            print(f"Error processing metric {metric}: {e}")

    try:
        metrics_to_average = ["NMI", "ARI", "ASW", "HOM", "COM"]
        metric_dfs = []
        for metric in metrics_to_average:
            if metric in metric_files:
                data = pd.read_csv(metric_files[metric], index_col=0)
                data_subset = data[dataset_names]
                metric_dfs.append(data_subset)

        # Compute composite score
        composite_data = sum(metric_dfs) / len(metrics_to_average)

        sns.heatmap(
            composite_data,
            cmap="coolwarm",
            cbar_kws={'label': "composite"},
            annot=True,  # Add annotations
            fmt=".2f",   # Format for annotations
            annot_kws={"size": 6 ,"weight": "bold" },  # Small font size for annotations
            ax=axes[-1]
        )
        ax.set_title(f"{dataset_type} Dataset Heatmap - Composite" , fontweight = "bold")
        ax.set_xlabel("Datasets" , fontweight = "bold")
        ax.set_ylabel("Methods" , fontweight = "bold")

        labels = ax.get_xticklabels()
        bold_font = fm.FontProperties(weight='bold')
        ax.set_xticklabels(labels, rotation=45, ha="right", fontproperties=bold_font)

        y_labels = ax.get_yticklabels()
        ax.set_yticklabels(y_labels, fontproperties=bold_font)

    except Exception as e:
        axes[-1].set_visible(False)  # Hide composite score plot if there's an error
        print(f"Error computing composite score: {e}")

 

    plt.tight_layout()  # Ensure spacing between plots
    plt.savefig(output_file, bbox_inches="tight", dpi = 200)
    plt.close()

def plot_individual_heatmaps(metric_files, dataset_names, dataset_type, output_dir):

    """
    Produce individual heatmaps for metrics of choice

    Parameters
    ----------
    metric_files : dictionary
        Mapping of metric names to output files 
    dataset_names: list
        Names of datasets to include in heatmap
    dataset_type: str
        Represents name of selected group of datasets (for example : simulated/real/DLPFC)
        Used only for naming
    output_file : str
        Path to save output image to
    """
    os.makedirs(output_dir, exist_ok=True)
    metrics = ["ARI", "NMI", "CHAOS", "PAS", "ASW", "HOM", "COM"]
    metric_dfs = []

    for metric, file_path in metric_files.items():
        try:
            data = pd.read_csv(file_path, index_col=0)
            data_subset = data[dataset_names]

            plt.figure(figsize=(8, 6))
            ax = sns.heatmap(
                data_subset,
                cmap="coolwarm",
                cbar_kws={'label': metric},
                annot=True,
                fmt=".2f",
                annot_kws={"size": 6, "weight": "bold"}
            )
            ax.set_title(f"{dataset_type} Dataset Heatmap - {metric}", fontweight="bold")
            ax.set_xlabel("Datasets", fontweight="bold")
            ax.set_ylabel("Methods", fontweight="bold")

            labels = ax.get_xticklabels()
            bold_font = fm.FontProperties(weight='bold')
            ax.set_xticklabels(labels, rotation=45, ha="right", fontproperties=bold_font)

            y_labels = ax.get_yticklabels()
            ax.set_yticklabels(y_labels, fontproperties=bold_font)

            output_path = os.path.join(output_dir, f"{dataset_type}_Heatmap_{metric}.png")
            plt.tight_layout()
            plt.savefig(output_path, bbox_inches="tight")
            plt.close()

            if metric in metrics:
                metric_dfs.append(data_subset)

        except Exception as e:
            print(f"Error processing metric {metric}: {e}")


    try:
        if metric_dfs:
            composite_data_accuracy = sum(metric_dfs[0:4]) / 4
            composite_data_accuracy.to_csv(os.path.join(output_dir, f"{dataset_type}_Heatmap_Composite_accuracy.csv"))
            plt.figure(figsize=(8, 6))
            ax = sns.heatmap(
                composite_data_accuracy,
                cmap="coolwarm",
                cbar_kws={'label': "composite - accuracy"},
                annot=True,
                fmt=".2f",
                annot_kws={"size": 6, "weight": "bold"}
            )
            ax.set_title(f"{dataset_type} Dataset Heatmap - Composite (accuracy)", fontweight="bold")
            ax.set_xlabel("Datasets", fontweight="bold")
            ax.set_ylabel("Methods", fontweight="bold")

            labels = ax.get_xticklabels()
            bold_font = fm.FontProperties(weight='bold')
            ax.set_xticklabels(labels, rotation=45, ha="right", fontproperties=bold_font)

            y_labels = ax.get_yticklabels()
            ax.set_yticklabels(y_labels, fontproperties=bold_font)

            output_path = os.path.join(output_dir, f"{dataset_type}_Heatmap_Composite_accuracy.png")
            plt.tight_layout()
            plt.savefig(output_path, bbox_inches="tight", dpi=200)
            plt.close()
    except Exception as e:
        print(f"Error computing composite accuracy score: {e}")

    try:
        if metric_dfs:
            composite_data_consistency = sum(metric_dfs[4:7]) / 3
            composite_data_consistency.to_csv(os.path.join(output_dir, f"{dataset_type}_Heatmap_Composite_consistency.csv"))

            plt.figure(figsize=(8, 6))
            ax = sns.heatmap(
                composite_data_consistency,
                cmap="coolwarm",
                cbar_kws={'label': "composite - consistency"},
                annot=True,
                fmt=".2f",
                annot_kws={"size": 6, "weight": "bold"}
            )
            ax.set_title(f"{dataset_type} Dataset Heatmap - Composite (consistency)", fontweight="bold")
            ax.set_xlabel("Datasets", fontweight="bold")
            ax.set_ylabel("Methods", fontweight="bold")

            labels = ax.get_xticklabels()
            bold_font = fm.FontProperties(weight='bold')
            ax.set_xticklabels(labels, rotation=45, ha="right", fontproperties=bold_font)

            y_labels = ax.get_yticklabels()
            ax.set_yticklabels(y_labels, fontproperties=bold_font)

            output_path = os.path.join(output_dir, f"{dataset_type}_Heatmap_Composite_consistency.png")
            plt.tight_layout()
            plt.savefig(output_path, bbox_inches="tight", dpi=200)
            plt.close()
    except Exception as e:
        print(f"Error computing composite consistency score: {e}")

    try:
        if metric_dfs:
            composite_data = (composite_data_accuracy + composite_data_consistency)/2
            composite_data.to_csv(os.path.join(output_dir, f"{dataset_type}_Heatmap_Composite.csv"))

            plt.figure(figsize=(8, 6))
            ax = sns.heatmap(
                composite_data,
                cmap="coolwarm",
                cbar_kws={'label': "composite"},
                annot=True,
                fmt=".2f",
                annot_kws={"size": 6, "weight": "bold"}
            )
            ax.set_title(f"{dataset_type} Dataset Heatmap - Composite", fontweight="bold")
            ax.set_xlabel("Datasets", fontweight="bold")
            ax.set_ylabel("Methods", fontweight="bold")

            labels = ax.get_xticklabels()
            bold_font = fm.FontProperties(weight='bold')
            ax.set_xticklabels(labels, rotation=45, ha="right", fontproperties=bold_font)

            y_labels = ax.get_yticklabels()
            ax.set_yticklabels(y_labels, fontproperties=bold_font)

            output_path = os.path.join(output_dir, f"{dataset_type}_Heatmap_Composite.png")
            plt.tight_layout()
            plt.savefig(output_path, bbox_inches="tight", dpi=200)
            plt.close()
    except Exception as e:
        print(f"Error computing composite score: {e}")


def make_embedding_plots(dataset_name,dataset_paths,pred_paths,save_dir, method_names = ["ground_truth","SCANIT", "CCST" , "DeepST" , "GraphST" , "PROST" , "SpaSRL" , "STAGATE","SpatialPCA" , "banksy" , "giotto" , "DR_SC" , "ISC_MEB" , "BayesSpace" , "PRECAST" , "BayesCafe" , "BASS"],stitch=True):
    """
    Produce embedding plots for select datasets and method outputs

    Parameters
    ----------
    dataset_name : str
        Name of dataset to create plot for
    dataset_paths : dict
        Dictionary of file paths needed for that dataset.
        Must contain at least 'st_path'.
        Can contain 'gnd_path' if ground truth is from a CSV/TSV file.
    pred_paths: dict of dicts
        Dictionary of output file paths for those datasets and names
        saved as a dataset x method dictionary
    save_dir : str
        Path to save output images to
    """
    
    method_names_dict = {
    "SCANIT": "SCANIT",
    "CCST": "CCST",
    "DeepST": "DeepST",
    "GraphST": "GraphST",
    "PROST": "PROST",
    "SpaSRL": "SpaSRL",
    "STAGATE": "STAGATE",
    "SpatialPCA": "SpatialPCA",
    "banksy": "Banksy",
    "giotto": "Giotto",
    "DR_SC": "DR.SC",
    "ISC_MEB": "ISC.MEB",
    "BayesSpace": "BayesSpace",
    "PRECAST": "PRECAST",
    "BayesCafe": "BayesCafe",
    "BASS": "BASS",
    "SpaceFlow":"SpaceFlow",
    "IRIS" : "IRIS"
    }

    dataset_names_dict = {
    "DLPFC151507": "DLPFC 151507",
    "DLPFC151508": "DLPFC 151508",
    "DLPFC151509": "DLPFC 151509",
    "DLPFC151510": "DLPFC 151510",
    "DLPFC151669": "DLPFC 151669",
    "DLPFC151670": "DLPFC 151670",
    "DLPFC151671": "DLPFC 151671",
    "DLPFC151672": "DLPFC 151672",
    "DLPFC151673": "DLPFC 151673",
    "DLPFC151674": "DLPFC 151674",
    "DLPFC151675": "DLPFC 151675",
    "DLPFC151676": "DLPFC 151676",
    "embryo9.5": "Embryo 9.5",
    "embryo14.5": "Embryo 14.5",
    "mouse_breast_cancer": "Mouse Breast Cancer",
    "MERFISH_brain0.04": "MERFISH Brain 0.04",
    "MERFISH_brain0.09": "MERFISH Brain 0.09",
    "MERFISH_brain0.14": "MERFISH Brain 0.14",
    "MERFISH_brain0.19": "MERFISH Brain 0.19",
    "MERFISH_brain0.24": "MERFISH Brain 0.24",
    "osmFISH": "osmFISH",
    "simulated_kidney_cancer410": "Kidney Cancer 410",
    "simulated_kidney_cancer411": "Kidney Cancer 411",
    "simulated_kidney_cancer506": "Kidney Cancer 506",
    "simulated_breast_cancerER+_CID4290": "Breast Cancer ER+ CID 4290",
    "simulated_breast_cancerTNBC_CID44971": "Breast Cancer TNBC CID 44971",
    "simulated_liver_cancerHCC-1L": "Liver Cancer HCC-1L",
    "simulated_liver_cancerHCC-2L": "Liver Cancer HCC-2L",
    "simulated_liver_cancerHCC-3L": "Liver Cancer HCC-3L",
    "simulated_liver_cancerHCC-4L": "Liver Cancer HCC-4L",
    "simulated_breast_atlas": "Breast Atlas",
    "simulated_intestineA1": "Intestine A1",
    "simulated_intestineA2": "Intestine A2",
    "simulated_chicken_heart": "Chicken Heart",
    "simulated_prostate_cancer" : "Prostate Cancer",
    "simulated_cerebellum" : "Cerebellum"
    }



    gnd, locs = import_dataset(dataset_name, dataset_paths[dataset_name], mode = "plot")

    for method_name in method_names:
        if method_name == "ground_truth":
            gnd = gnd.squeeze()
            adata = an.AnnData(obs=pd.DataFrame({"cluster": gnd.astype(str)}))
            adata.obsm["spatial"] = np.array(locs)
        else:
            pred = load_prediction(dataset_name,method_name, pred_paths[dataset_name][method_name])
            if pred is None:
                print(f"output for {dataset_name} {method_name} not present")
                continue
            pred = pred[~pred.index.duplicated(keep='first')]
            intersect_idx = gnd.index.intersection(pred.index)
            pred_filtered = pred.loc[intersect_idx]

            # Handle both DataFrame and Series cases safely
            if isinstance(pred_filtered, pd.DataFrame):
                if pred_filtered.shape[1] > 1:
                    # If multiple columns (e.g. banksy), take the first or warn
                    print(f"Warning: {method_name} prediction for {dataset_name} has multiple columns; using the first one.")
                pred_final = pred_filtered.iloc[:, 0]
            else:
                pred_final = pred_filtered

            pred_final = pred_final.squeeze()  # Ensure Series, not DataFrame
            pred_final.index = intersect_idx   # Ensure correct indexing
            # pred_values = pred_filtered.iloc[:, 0].values
            spatial_coords = locs.loc[intersect_idx].values
            adata = an.AnnData(obs=pd.DataFrame({"cluster": pred_final.astype(str)}))
            adata.obsm['spatial'] = np.array(spatial_coords)
                

        # Ensure the directory exists
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
            print(f"Created directory: {save_dir}")
    
        # Generate and save the plot
        plt.rcParams['font.weight'] = 'bold'       # Bold for all text
        plt.rcParams['axes.titleweight'] = 'bold'  # Bold for titles
        plt.rcParams['axes.labelweight'] = 'bold'
        fig = sc.pl.embedding(
            adata,
            basis="spatial",  # Use the 'spatial' embedding
            color="cluster",  # Color points by the 'cluster' column
            title=f"{dataset_names_dict[dataset_name]} - {"Ground Truth" if method_name == "ground_truth" else method_names_dict[method_name]}",
            size=100,
            alpha=1,
            cmap="coolwarm",
            show = True,
            return_fig=True
            # save=f"{dataset_name}_{method_name}.png"  # File will be saved in the default directory by Scanpy
        )

        
        # Move the file to the desired directory
        current_dir = os.getcwd()
        plot_path = os.path.join(current_dir, f"figures/spatial{dataset_name}.png")
        if not os.path.exists(os.path.join(save_dir, f"{dataset_name}")):
            os.makedirs(os.path.join(save_dir, f"{dataset_name}"))
        new_path = os.path.join(save_dir, f"{dataset_name}/{method_name}.png")
        fig.savefig(new_path, dpi=200, bbox_inches="tight", facecolor="white")
        print(f"figure saved to {new_path}")
        plt.close(fig)

    
        # if os.path.exists(plot_path):
        #     os.rename(plot_path, new_path)
        #     print(f"Moved plot to: {new_path}")
        # else:
        #     print(f"Plot file not found at: {plot_path}")
        
    
    
# Helper functions to stitch embedding plots together
def get_all_files(folder_path):
    file_names = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
    png_files = [f for f in file_names if f.lower().endswith('.png')]
    png_files.sort()

    ground_truth = 'ground_truth.png'  # Replace with the actual name or pattern for the ground truth file
    if ground_truth in png_files:
        png_files.remove(ground_truth)
        png_files.insert(0, ground_truth)
    return png_files

def save_grid_image(grid_rows, grid_cols, folder_path, dataset_name, save_path, width, height):
    png_files = get_all_files(folder_path)
    image_paths = [os.path.join(folder_path, p) for p in png_files]
    
    images = [Image.open(img) for img in image_paths]
    fig, axes = plt.subplots(grid_rows, grid_cols, figsize=(width, height), facecolor="white")
    fig.subplots_adjust(wspace=0.0, hspace=0.0)  # Reduced spacing between rows and columns
    axes = axes.flatten()
    for ax, img in zip(axes, images + [None] * (grid_rows * grid_cols - len(images))):
        if img is not None:
            ax.imshow(img)
        ax.axis("off")  # Remove axes

    output_filename = dataset_name + "_grid_image.png"
    plt.savefig(os.path.join(save_path, output_filename), dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Grid saved as {output_filename} with 200 DPI.")

# ----------------------- Composite Score Computation (over data subsets) -----------------------


def return_custom_score(metric_files, dataset_names):
    metric_dfs = []

    for metric, file_path in metric_files.items():
        try:
            data = pd.read_csv(file_path, index_col=0)
            data_subset = data[dataset_names]
            metric_dfs.append(data_subset)

        except Exception as e:
            print(f"Error processing metric {metric}: {e}")


    try:
        if metric_dfs:
            composite_data_accuracy = sum(metric_dfs[0:4]) / 4
            
    except Exception as e:
        print(f"Error computing composite accuracy score: {e}")

    try:
        if metric_dfs:
            composite_data_consistency = sum(metric_dfs[4:7]) / 3
            

    except Exception as e:
        print(f"Error computing composite consistency score: {e}")

    try:
        if metric_dfs:
            composite_data = (composite_data_accuracy + composite_data_consistency)/2
            
    except Exception as e:
        print(f"Error computing composite score: {e}")

    composite_data = np.sum(composite_data, axis = 1)/len(dataset_names)
    top_5_rows = (
        composite_data             # Convert to Series with MultiIndex (row, column)
        .sort_values(ascending=False)  # Sort in descending order
        .head(5)               # Get top 5
        .index.get_level_values(0)  # Extract row names
        .tolist()              # Convert to list
    )


    return(composite_data, top_5_rows)


# ----------------------- Useful Variables (comment out to use) -----------------------
# dataset_names = ["DLPFC151507" , "DLPFC151508" , "DLPFC151509","DLPFC151510","DLPFC151669","DLPFC151670","DLPFC151671","DLPFC151672","DLPFC151673","DLPFC151674","DLPFC151675","DLPFC151676","embryo9.5" , "embryo14.5" , "mouse_breast_cancer", "MERFISH_brain0.04", "MERFISH_brain0.09", "MERFISH_brain0.14", "MERFISH_brain0.19", "MERFISH_brain0.24","osmFISH" , "simulated_kidney_cancer410" , "simulated_kidney_cancer411" , "simulated_kidney_cancer506" , "simulated_breast_cancerER+_CID4290"  , "simulated_breast_cancerTNBC_CID44971" , "simulated_liver_cancerHCC-1L", "simulated_liver_cancerHCC-2L", "simulated_liver_cancerHCC-3L", "simulated_liver_cancerHCC-4L" , "simulated_breast_atlas" , "simulated_intestineA1", "simulated_intestineA2","simulated_chicken_heart", "simulated_prostate_cancer", "simulated_cerebellum"]

# DLPFC_dataset_names = ["DLPFC151507" , "DLPFC151508" , "DLPFC151509","DLPFC151510","DLPFC151669","DLPFC151670","DLPFC151671","DLPFC151672","DLPFC151673","DLPFC151674","DLPFC151675","DLPFC151676"]

# real_dataset_names = ["embryo9.5" , "embryo14.5" , "mouse_breast_cancer", "MERFISH_brain0.04", "MERFISH_brain0.09", "MERFISH_brain0.14", "MERFISH_brain0.19", "MERFISH_brain0.24","osmFISH"]

# simulated_dataset_names = ["simulated_kidney_cancer410" , "simulated_kidney_cancer411" , "simulated_kidney_cancer506" , "simulated_breast_cancerER+_CID4290"  , "simulated_breast_cancerTNBC_CID44971" , "simulated_liver_cancerHCC-1L", "simulated_liver_cancerHCC-2L", "simulated_liver_cancerHCC-3L", "simulated_liver_cancerHCC-4L" , "simulated_breast_atlas" , "simulated_intestineA1", "simulated_intestineA2","simulated_chicken_heart", "simulated_prostate_cancer", "simulated_cerebellum"]

# method_names = ["SCANIT", "CCST" , "DeepST" , "GraphST" , "PROST" , "SpaSRL" , "STAGATE","SpatialPCA" , "banksy" , "giotto" , "DR_SC" , "ISC_MEB" , "BayesSpace" , "PRECAST" , "BayesCafe" , "BASS","SpaceFlow", "IRIS"]



