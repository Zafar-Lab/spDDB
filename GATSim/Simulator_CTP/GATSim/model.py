import tensorflow.compat.v1 as tf
import tensorflow as tf2
import numpy as np

tf.compat.v1.enable_eager_execution()

class GATE():

    def __init__(self, hidden_dims, alpha=0.8, nonlinear=True, weight_decay=0.0001):
        
        self.n_layers = len(hidden_dims) - 1
        self.alpha = alpha
        self.W, self.v, self.prune_v, self.W_logvar, self.W_mean, self.Ws_att_logvar, self.Ws_att_mean = self.define_weights(hidden_dims)
        self.C = {}
        self.prune_C = {}
        self.nonlinear = nonlinear
        self.weight_decay = weight_decay

    def reconstruction(self, X, X_alpha):
                
        dist = tf.distributions.Dirichlet(X_alpha)
        neg_loglikelihood = tf.reduce_sum(-1*dist.log_prob(X))
        
        return neg_loglikelihood

    def reparameterize(self, mean, variance, scope=None):

        with tf.variable_scope(scope, 'sample_gaussian'):
            sample = tf.random_normal(tf.shape(mean), mean, tf.exp(variance))
            sample.set_shape(mean.get_shape())
            return sample
            
    def dense(self, x, inp_dim, out_dim, name = 'dense'):

        with tf.variable_scope(name, reuse=None):
            weights = tf.get_variable("weights", shape=[inp_dim, out_dim],
                                      initializer =
                                      tf2.initializers.GlorotUniform()) 
            
            bias = tf.get_variable("bias", shape=[out_dim], initializer = tf.constant_initializer(0.0))
            
            out = tf.add(tf.matmul(x, weights), bias, name='matmul')
            return out

    def last_layer(self, A, H, layer):

        H_mu = tf.matmul(H, self.W_mean)
        H_logvar = tf.matmul(H, self.W_logvar)

        self.C[2] = self.graph_attention_layer(A, H_mu, self.Ws_att_mean, layer)
        H_mu = tf.sparse_tensor_dense_matmul(self.C[2], H_mu)
        H_logvar = tf.sparse_tensor_dense_matmul(self.graph_attention_layer(A, H_mu,  self.Ws_att_logvar, layer + 1), H_logvar)

        return H_mu, H_logvar
    
    def __call__(self, A, prune_A, X, A_dense_reshaped):
        
        # Encoder
        H = X

        for layer in range(self.n_layers): # AJ


            H = self.__encoder(A, prune_A, H, layer)
            if self.nonlinear:

                if layer != self.n_layers-1: # Removed as adding GAT layer.
                    H = tf.nn.elu(H)

        self.H = H
        
        for layer in range(self.n_layers - 1, -1, -1):
            H = self.__decoder(H, layer)
            if self.nonlinear:
                if layer != 0:
                    H = tf.nn.elu(H)
        X_ = H

        # Reconstruction loss for Adjacency matrix
        A_pred = tf.nn.sigmoid(tf.matmul(self.H, tf.transpose(self.H)))    
        self.A_pred_loss = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(tf.reshape(A_dense_reshaped, [-1]), tf.reshape(A_pred, [-1])))

        # Reconstruction loss for Cell type proportion matrix
        features_loss = tf.sqrt(tf.reduce_sum(tf.reduce_sum(tf.pow(X - X_, 2))))

        self.kl_loss = tf.constant(0, dtype = tf.float32)

        for layer in range(self.n_layers):
            weight_decay_loss = 0
            weight_decay_loss += tf.multiply(tf.nn.l2_loss(self.W[layer]), self.weight_decay, name = 'weight_loss')
        
        self.loss = features_loss  + weight_decay_loss + (self.kl_loss + self.A_pred_loss)
        
        if self.alpha == 0:
            self.Att_l = self.C
        else:
            self.Att_l = {'C': self.C, 'prune_C': self.prune_C}
        
        return self.loss, self.kl_loss, self.A_pred_loss, self.H, self.Att_l, X_

    def __encoder(self, A, prune_A, H, layer):
        H = tf.matmul(H, self.W[layer])
        
        if layer == self.n_layers-1: 
            return H
        
        self.C[layer] = self.graph_attention_layer(A, H, self.v[layer], layer)
        
        if self.alpha == 0:
            return tf.sparse_tensor_dense_matmul(self.C[layer], H)
        
        else:
            self.prune_C[layer] = self.graph_attention_layer(prune_A, H, self.prune_v[layer], layer)
            
            return (1-self.alpha)*tf.sparse_tensor_dense_matmul(self.C[layer], H) 
            + self.alpha*tf.sparse_tensor_dense_matmul(self.prune_C[layer], H)


    def __decoder(self, H, layer):
        H = tf.matmul(H, self.W[layer], transpose_b=True)
        if layer == 0:
            return H
        if self.alpha == 0:
            return tf.sparse_tensor_dense_matmul(self.C[layer-1], H)
        else:
            return (1-self.alpha)*tf.sparse_tensor_dense_matmul(self.C[layer-1], H) 
            + self.alpha*tf.sparse_tensor_dense_matmul(self.prune_C[layer-1], H)
            
    def define_weights(self, hidden_dims):
        W = {}
        for i in range(self.n_layers): #AJ -1 
            W[i] = tf.get_variable("W%s" % i, shape=(hidden_dims[i], hidden_dims[i+1]))

        W_logvar = tf.get_variable("W_logvar", shape=(10, 10))
        W_mean = tf.get_variable("W_mean", shape=(10, 10))
        
        Ws_att = {}

        for i in range(self.n_layers-1): # AJ -1
            v = {}
            v[0] = tf.get_variable("v%s_0" % i, shape=(hidden_dims[i+1], 1))
            v[1] = tf.get_variable("v%s_1" % i, shape=(hidden_dims[i+1], 1))

            Ws_att[i] = v
            
        
        Ws_att_logvar = {}
        Ws_att_logvar[0] = tf.get_variable("v%s_0logvar" % i, shape=(10, 1))
        Ws_att_logvar[1] = tf.get_variable("v%s_1logvar" % i, shape=(10, 1))

        Ws_att_mean = {}
        Ws_att_mean[0] = tf.get_variable("v%s_0mean" % i, shape=(10, 1))
        Ws_att_mean[1] = tf.get_variable("v%s_1mean" % i, shape=(10, 1))
        
        
        prune_Ws_att = {}
        
        if self.alpha == 0:
            return W, Ws_att, prune_Ws_att, W_logvar, W_mean, Ws_att_logvar, Ws_att_mean
            #return W, Ws_att, prune_Ws_att
        
        for i in range(self.n_layers-1):
            prune_v = {}
            prune_v[0] = tf.get_variable("prune_v%s_0" % i, shape=(hidden_dims[i+1], 1))
            prune_v[1] = tf.get_variable("prune_v%s_1" % i, shape=(hidden_dims[i+1], 1))

            prune_Ws_att[i] = prune_v

        return W, Ws_att, prune_Ws_att, W_logvar, W_mean, Ws_att_logvar, Ws_att_mean

    def graph_attention_layer(self, A, M, v, layer):

        with tf.variable_scope("layer_%s"% layer):
            f1 = tf.matmul(M, v[0])
            f1 = A * f1

            print ("f1 shape is:", f1.shape)
            
            f2 = tf.matmul(M, v[1])
            f2 = A * tf.transpose(f2, [1, 0])
            logits = tf.sparse_add(f1, f2)

            unnormalized_attentions = tf.SparseTensor(indices=logits.indices,
                                         values=tf.nn.sigmoid(logits.values),
                                         dense_shape=logits.dense_shape)
            attentions = tf.sparse_softmax(unnormalized_attentions)

            attentions = tf.SparseTensor(indices=attentions.indices,
                                         values=attentions.values,
                                         dense_shape=attentions.dense_shape)

            return attentions