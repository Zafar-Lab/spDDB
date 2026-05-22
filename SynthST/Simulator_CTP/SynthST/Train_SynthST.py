import numpy as np
import scipy.sparse as sp
import tensorflow.compat.v1 as tf
import pandas as pd
import scanpy as sc
from SynthST.SynthST import SynthST 
from sklearn.preprocessing import normalize

def train_SynthST(adata_list, method_list, hidden_dims=[512, 30], alpha=0, n_epochs=1000, lr=0.0001, key_added='SynthST',
                gradient_clipping=5, nonlinear=True, weight_decay=0.0001,verbose=True, 
                random_seed=2020, pre_labels=None, pre_resolution=0.2,
                save_attention=False, save_loss=False, save_reconstrction=False):

    tf.compat.v1.disable_eager_execution()
    tf.reset_default_graph()
    np.random.seed(random_seed)
    tf.set_random_seed(random_seed)

    adata = adata_list[0]
    if 'highly_variable' in adata.var.columns:
        adata_Vars =  adata[:, adata.var['highly_variable']]
    else:
        adata_Vars = adata

    X = pd.DataFrame(adata_Vars.X[:, ], index=adata_Vars.obs.index, columns=adata_Vars.var.index)
    cells = np.array(X.index)
    
    trainer = SynthST(hidden_dims=[X.shape[1]] + hidden_dims, alpha=alpha, 
                    n_epochs=n_epochs, lr=lr, gradient_clipping=gradient_clipping, 
                    nonlinear=nonlinear,weight_decay=weight_decay, verbose=verbose, 
                    random_seed=random_seed)

    if alpha == 0:
        
        X_list_testing = []
        for i in range(len(adata_list)):

            adata = adata_list[i]
            if 'highly_variable' in adata.var.columns:
                adata_Vars =  adata[:, adata.var['highly_variable']]
            else:
                adata_Vars = adata
    
            X = pd.DataFrame(adata_Vars.X[:, ], index=adata_Vars.obs.index, columns=adata_Vars.var.index) #.toarray()
            X_list_testing += [X]


        trainer(adata_list) ## G_tf , G_tf, X earlier
        index = 0

        arrays_stacked = np.stack(X_list_testing)
        average_array = pd.DataFrame(np.mean(arrays_stacked, axis=0), index = adata_Vars.obs.index, columns = adata_Vars.var.index)
        X_list_testing += [average_array]

        median_array = pd.DataFrame(np.mean(arrays_stacked, axis=0), index = adata_Vars.obs.index, columns = adata_Vars.var.index)
        X_list_testing += [median_array]
        
        
        for X in X_list_testing:
            print ("Inference for :", method_list[index])
            
            #embeddings, attentions, loss, ReX, X_alpha = trainer.infer(X)
            embeddings, attentions, loss, ReX = trainer.infer(X)
            
            cell_reps = pd.DataFrame(embeddings)
            cell_reps.index = cells
        
            adata.obsm[method_list[index] + "_embedding"] = cell_reps.loc[adata.obs_names, ].values
            
            if save_attention:
                adata.uns['SynthST_attention'] = attentions
                
            if save_loss:
                adata.uns['SynthST_loss'] = loss
                
            if save_reconstrction:
                           
                ReX = pd.DataFrame(ReX, index=X.index, columns=X.columns)
                ReX[ReX<0] = 0
               
                # We will only store Norm so that the final CTP matrix is always normalized
                #adata.obsm[method_list[index] + "_SynthST_ReX"] = ReX.values

                ReX = ReX.div(ReX.sum(axis = 1), axis = 0)
                print (np.sort(ReX.sum(axis = 1)))
                adata.obsm[method_list[index] + "_SynthST_ReX_Norm"] = ReX #.values
                
            index += 1
   
    return adata


def prune_spatial_Net(Graph_df, label):
    print('------Pruning the graph...')
    print('%d edges before pruning.' %Graph_df.shape[0])
    pro_labels_dict = dict(zip(list(label.index), label))
    Graph_df['Cell1_label'] = Graph_df['Cell1'].map(pro_labels_dict)
    Graph_df['Cell2_label'] = Graph_df['Cell2'].map(pro_labels_dict)
    Graph_df = Graph_df.loc[Graph_df['Cell1_label']==Graph_df['Cell2_label'],]
    print('%d edges after pruning.' %Graph_df.shape[0])
    return Graph_df


def recovery_Imputed_Count(adata, size_factor):
    assert('ReX' in adata.uns)
    temp_df = adata.uns['ReX'].copy()
    sf = size_factor.loc[temp_df.index]
    temp_df = np.expm1(temp_df)
    temp_df = (temp_df.T * sf).T
    adata.uns['ReX_Count'] = temp_df
    return adata
