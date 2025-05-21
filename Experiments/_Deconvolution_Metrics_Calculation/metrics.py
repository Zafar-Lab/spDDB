import pandas as pd
import numpy as np
import scanpy as sc
import scipy
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import entropy
from scipy.spatial.distance import jensenshannon
import seaborn as sns
import os
import pickle
from matplotlib import rcParams
from scipy.spatial import distance
from sklearn.metrics import precision_recall_curve, auc

def get_weight_matrix(coords,l = 1,co = 0, eps = 0.00000001):
    
    distances = squareform(pdist(coords, metric='euclidean'))
    weight_matrix = np.exp((-0.5)*distances**2/l**2)
    sums = np.sum(weight_matrix)
    weight_matrix = len(coords)*weight_matrix/sums
    #weight_matrix[weight_matrix < co] = 0
   # weight_matrix = np.eye((len(coords)))
    return weight_matrix



def get_morans_R(x,y,coords, l = 1, co = 0, eps = 0.00000001):
    w = get_weight_matrix(coords, l = l, co = co)
    x = np.array(x)
    y = np.array(y)
    
    #print ("x", x, "y", y)
    
    if(np.all(y==0)):
        return 0
    x_bar = np.mean(x)
    y_bar = np.mean(y)
    x = x - x_bar
    y = y - y_bar
    l_x = np.sqrt(np.sum(x*x))
    l_y = np.sqrt(np.sum(y*y))

    x = x.reshape(-1,1)
    y = y.reshape(-1,1)
    return (np.sum(w*(x@y.T))/(l_x*l_y + eps)+1)/2


def get_cosine_sim(x,y, eps = 0.00000001):
    #w = get_weight_matrix(coords, l = l, co = co)
    x = np.array(x)
    y = np.array(y)
    if(np.all(y==0)):
        return 0
    return (2-distance.cosine(x,y))/2
    n1 = (w*x*y).sum(axis=1)
    n2 = (w*x*x).sum(axis=1)
    n3 = (w*y*y).sum(axis=1)
    return (np.mean(n1/(np.sqrt(n2*n3) + eps))+1)/2
    
def get_spatial_pearson(x,y,coords, l = 1, co = 0, eps = 0.00000001):

    w = get_weight_matrix(coords, l = l, co = co)
    x = np.array(x)
    y = np.array(y)
    
    if(np.all(y==0)):
        return 0
 
    n1 = w.sum(axis=1)
    n2 = (w*x*y).sum(axis=1)
    n3 = (w*x).sum(axis=1)
    n4 = (w*y).sum(axis=1)
    num = n1*n2 - n3*n4
    n2 = (w*x*x).sum(axis=1)
    n3 = (w*y*y).sum(axis=1)
    n4 = (w*x).sum(axis=1)
    n5 = (w*y).sum(axis=1)
    n4 = n4*n4
    n5 = n5*n5
    
    denom = (n1*n2 - n4)*(n1*n3-n5)
    
    if (np.any(denom < 0)):
        print ("Debug", denom[denom < 0], "\n")
        denom[denom < 0] = eps

    """
    if (np.any(np.isnan(denom))):
        print ("Nan in denom")
    elif (np.any(np.isinf(denom))):
        print ("Infinity in denom")
    """    
    den = np.sqrt(denom) + eps
    return (np.mean(num/den)+1)/2


def get_spearman(x,y,coords, l = 1, co = 0, eps = 0.00000001):
    w = get_weight_matrix(coords, l = l, co = co)
    x = np.array(x)
    y = np.array(y)
  
    if(np.all(y==0)):
        return 0
    x = np.argsort(x)
    y = np.argsort(y)
    
    n1 = w.sum(axis=1)
    n2 = (w*x*y).sum(axis=1)
  
    n3 = (w*x).sum(axis=1)
    n4 = (w*y).sum(axis=1)
    num = n1*n2 - n3*n4
    n2 = (w*x*x).sum(axis=1)
    n3 = (w*y*y).sum(axis=1)
    n4 = (w*x).sum(axis=1)
    n5 = (w*y).sum(axis=1)
    n4 = n4*n4
    n5 = n5*n5

    den = np.sqrt((n1*n2 - n4)*(n1*n3-n5)) + eps
    return (np.mean(num/den)+1)/2

"""
def get_rmse(preds,gt, eps = 0.00000001):
    rmse = np.sqrt((np.array(preds)-np.array(gt))**2)
    return 1-np.mean(rmse)


def get_rmse(preds, gt, eps = 0.00000001):

    rmse_list = []
    preds = np.array(preds)
    gt = np.array(gt)
    
    for i in range(len(preds)):
        rmse = np.sqrt(np.mean((gt[i, :] - preds[i, :]) ** 2))
        rmse_list.append(rmse)
        
    return 1-(sum(rmse_list)/len(preds)) 
"""

def get_rmse(preds, gt, eps = 0.00000001):

    rmse_list = []
    preds = np.array(preds)
    gt = np.array(gt)
    
    for i in range(len(preds)):
        rmse = np.sum((gt[i, :] - preds[i, :]) ** 2)
        rmse_list.append(rmse)
        
    return (1 - np.sqrt(sum(rmse_list)/len(preds))) 


def get_jsd(preds,gt, eps = 0.00000001):
    jsd_list = []
    preds = np.array(preds)
    gt = np.array(gt)
    for i in range(len(preds)):
        jsd = jensenshannon(preds[i,:], gt[i,:]) # To compute row wise.
        if(jsd>100 or np.isnan(jsd)):
            jsd_list.append(1)
            continue
        jsd_list.append(jsd)
    return 1-sum(jsd_list)/len(preds) 
    
def get_pearson(x,y, eps = 0.00000001):
    x = np.array(x)
    y = np.array(y)
    
    if (x.std() == 0):
        print ("Ground truth's standard deviation is 0")
    elif (y.std() == 0):
        print ("Predictions' standard deviation is 0")
        
    if(np.all(y==0)):
        return 0
    return (np.corrcoef(x,y)[0, 1]+1)/2


def compute_Lee_stats(x, y, coords, l = 1, co = 0, eps = 0.00000001):
    
    w = get_weight_matrix(coords, l = l, co = co)
    x = np.array(x)
    y = np.array(y)
    
    if(np.all(y == 0)):
        return 0

    x_bar = np.mean(x)
    y_bar = np.mean(y)
    x = x - x_bar
    y = y - y_bar
    
    # Constant term...

    n = x.shape[0]
    w_sum = np.sum(w, axis = 1)
    const_term = n/np.sum(w_sum * w_sum) # Element wise sqaure.

    # Numerator...  
    
    #num = np.sum(np.sum((w @ x).T, axis = 1) * np.sum((w @ y).T, axis = 1))
    num = np.sum(w@x * w@y)
    
    # Denominator...
    
    l_x = np.sqrt(np.sum(x*x))
    l_y = np.sqrt(np.sum(y*y))

    denom = (l_x * l_y) + eps

    L_stat = const_term * (num/denom)
    return (L_stat+1)/2


    
def compute_ssim(x, y, eps = 0.00000001):
    
    C1 = 0.01
    C2 = 0.03

    #print ("x", x, " y", y)
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    sigma_x = np.std(x)
    sigma_y = np.std(y)
    
    cov = np.cov(x, y)[0, 1]

    #print ("cov", cov)
    
    Num = (2 * mean_x * mean_y + C1**2) * (2 * cov + C2**2)
    Denom = (mean_x**2 + mean_y**2 + C1**2) * (sigma_x**2 + sigma_y**2 + C2**2) + eps

    #print ("Num", Num, " Demon", Denom)
    
    SSIM = Num/(Denom + eps)
    #print (SSIM)
    return (SSIM+1)/2
    
def compute_geary(x , y ,coords, l = 1, co = 0, eps = 0.00000001):
  weights = get_weight_matrix(coords, l = l, co = co)
  x = np.array(x)
  y = np.array(y)
  weights = np.array(weights)
  S = np.sum(weights)
  x_bar = np.mean(x)
  y_bar = np.mean(y)
  denominator = 2*S*np.sqrt(np.sum((x-x_bar)**2))*np.sqrt(np.sum((y-y_bar)**2)) + eps
  numerator = 0
  for i in range(len(x)):
      numerator = numerator + np.sum(weights[i,:]*((x[i]-y)**2))
  numerator = (len(x)-1)*numerator
  
  return np.exp(-np.log(2)*(numerator/denominator))


def compute_AUPR(x, y ,coords, l = 1, co = 0, eps = 0.00000001):

    # x: gt, y:pred
    x = np.array(x)
    y = np.array(y)
    aupr_list = []
    
    thresholds = np.arange(0.01, 1.0, 0.1)

    print ("thresholds", thresholds)
    for thresh in thresholds:
    
        ground_truth = (x > thresh).astype(int)
        precision, recall, thresholds = precision_recall_curve(ground_truth, y)        
        # Compute AUPR (Area Under the Precision-Recall Curve)
        aupr = auc(recall, precision)
        aupr_list += [aupr]

    print (aupr_list)





























