import tensorflow.compat.v1 as tf
import scipy.sparse as sp
import numpy as np
import pandas as pd
from tqdm import tqdm
from GATSim.model import GATE
import random
import sys

class GATSim():

    def __init__(self, hidden_dims, alpha, n_epochs=500, lr=0.001, 
                 gradient_clipping=5, nonlinear=True, weight_decay=0.0001, 
                 verbose=True, random_seed=2020):
        
        np.random.seed(random_seed)
        tf.set_random_seed(random_seed)
        self.loss_list = []
        self.epoch_loss = [] # To take average of loss in each epoch
        self.lr = lr
        self.n_epochs = n_epochs
        self.gradient_clipping = gradient_clipping
        self.build_placeholders()
        self.verbose = verbose
        self.alpha = alpha
        self.gate = GATE(hidden_dims, alpha, nonlinear, weight_decay)
        self.loss, self.kl_loss, self.A_pred_loss, self.H, self.C, self.ReX = self.gate(self.A, self.prune_A, self.X, self.A_dense_reshaped)
        self.optimize(self.loss)
        self.build_session()
        tf.compat.v1.disable_eager_execution()
    
    def build_placeholders(self):
        self.A = tf.sparse_placeholder(dtype=tf.float32)
        self.prune_A = tf.sparse_placeholder(dtype=tf.float32)
        self.X = tf.placeholder(dtype=tf.float32)
        self.A_dense_reshaped = tf.placeholder(dtype=tf.float32)

    def build_session(self, gpu= True):
        config = tf.ConfigProto()
        config.gpu_options.allow_growth = True
        if gpu == False:
            config.intra_op_parallelism_threads = 0
            config.inter_op_parallelism_threads = 0
        self.session = tf.Session(config=config)
        self.session.run([tf.global_variables_initializer(), tf.local_variables_initializer()])

    def optimize(self, loss):
        optimizer = tf.train.AdamOptimizer(learning_rate=self.lr)
        gradients, variables = zip(*optimizer.compute_gradients(loss))
        gradients, _ = tf.clip_by_global_norm(gradients, self.gradient_clipping)
        self.train_op = optimizer.apply_gradients(zip(gradients, variables))

    def prepare_graph_data(self, adj):

        num_nodes = adj.shape[0]
        adj = adj + sp.eye(num_nodes)# self-loop
        
        #data =  adj.tocoo().data
        #adj[adj > 0.0] = 1.0
        
        if not sp.isspmatrix_coo(adj):
            adj = adj.tocoo()
        adj = adj.astype(np.float32)
        indices = np.vstack((adj.col, adj.row)).transpose()
        return (indices, adj.data, adj.shape)
    
    def __call__(self, adata_list):

        X_list = []
        no_of_samples = len(adata_list)
        
        #for adata in adata_list:
        for i in range(len(adata_list)):

            adata = adata_list[i]
            if 'highly_variable' in adata.var.columns:
                adata_Vars =  adata[:, adata.var['highly_variable']]
            else:
                adata_Vars = adata
    
            X = pd.DataFrame(adata_Vars.X[:, ], index=adata_Vars.obs.index, columns=adata_Vars.var.index) #.toarray()
            X_list += [X]
            
            if (i == (no_of_samples - 1)):
                if self.verbose:
                    print('Size of Input: ', adata_Vars.shape)
        
                cells = np.array(X.index)
                cells_id_tran = dict(zip(cells, range(cells.shape[0])))
        
                if 'Spatial_Net' not in adata.uns.keys():
                    raise ValueError("Spatial_Net is not existed! Run Cal_Spatial_Net first!")
        
                Spatial_Net = adata.uns['Spatial_Net']
        
                G_df = Spatial_Net.copy()
                G_df['Cell1'] = G_df['Cell1'].map(cells_id_tran)
                G_df['Cell2'] = G_df['Cell2'].map(cells_id_tran)
        
                G = sp.coo_matrix((np.ones(G_df.shape[0]), (G_df['Cell1'], G_df['Cell2'])), shape=(adata.n_obs, adata.n_obs))
                G_tf = self.prepare_graph_data(G)
                
                G_tf_dense = G.toarray()
    
                self.test_X = X
                self.test_G_tf = G_tf
                
        for epoch in tqdm(range(self.n_epochs), file = sys.stdout):
            
            np.random.shuffle(X_list)
            #print ("Xlist shuffling ", X_list)
            count = 0
            flag = False
            self.epoch_loss = []
            
            for X in X_list:
                
                #print (G_tf_dense) #AJ2024                  
                count += 1
                if (count == len(X_list)):
                    flag = True
                self.run_epoch(epoch, G_tf, G_tf, X, G_tf_dense, flag)
            
    def run_epoch(self, epoch, A, prune_A, X, A_dense_reshaped, flag):

        loss, _, kl_loss, A_pred_loss = self.session.run([self.loss, self.train_op, self.kl_loss, self.A_pred_loss],
                                 feed_dict = {self.A: A,
                                            self.prune_A: prune_A,
                                            self.X: X,
                                            self.A_dense_reshaped: A_dense_reshaped})                                          
        self.loss_list.append(loss)  
        self.epoch_loss.append(loss)
        #if self.verbose:
        if flag:
            print("Epoch: %s, Loss: %.4f" % (epoch, np.mean(self.epoch_loss)))
            
        return loss
        
    def infer(self, test_X): # A, prune_A, X commented

        H, C, ReX = self.session.run([self.H, self.C, self.ReX], # self.X_alpha
                           feed_dict = {self.A: self.test_G_tf,
                                      self.prune_A: self.test_G_tf,
                                      self.X:test_X}) # self.test_X

        return H, self.Conbine_Atten_l(C), self.loss_list, ReX #, X_alpha


    def Conbine_Atten_l(self, input):
        if self.alpha == 0:
            return [sp.coo_matrix((input[layer][1], (input[layer][0][:, 0], input[layer][0][:, 1])), shape=(input[layer][2][0], input[layer][2][1])) for layer in input]
        else:
            Att_C = [sp.coo_matrix((input['C'][layer][1], (input['C'][layer][0][:, 0], input['C'][layer][0][:, 1])), shape=(input['C'][layer][2][0], input['C'][layer][2][1])) for layer in input['C']]
            Att_pruneC = [sp.coo_matrix((input['prune_C'][layer][1], (input['prune_C'][layer][0][:, 0], input['prune_C'][layer][0][:, 1])), shape=(input['prune_C'][layer][2][0], input['prune_C'][layer][2][1])) for layer in input['prune_C']]
            return [self.alpha*Att_pruneC[layer] + (1-self.alpha)*Att_C[layer] for layer in input['C']]
