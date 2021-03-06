# -*- coding: utf-8 -*-
"""
Created on Sat Feb  9 18:01:40 2019

@author: ly
"""


 
#%%[1]input kmer.csv


def usage ():
    print("usage:")
    print("  Python encoder.py -i inputfile -n num_hidden -e iteration_number -l learning rate")
    print("    -i inputfile should be in csv format")
# In[1]:
import sys
import getopt

opts, args = getopt.getopt(sys.argv[1:], "hi:n:e:l:")
input_file=""

for op, value in opts:
    if op == "-i":
        input_file = value
    elif op == "-n":
        n_hidden=int(value)
    elif op == "-h":
        usage()
        sys.exit()
    elif op=="-e":
        training_epochs=int(value)
    elif op=="-l":
        learning_rate=float(value)
        
 
import numpy as np
 
import sklearn.preprocessing as prep
 
import tensorflow as tf
 
import csv

#input_file="Contigs_Sharon.txt_4mer.csv"

csv_reader=csv.reader(open(input_file))

title=[]
X=[]

for row in csv_reader:
    title.append(row[0])
    X.append(row[1:])
 
#%%[2]。
 

def xavier_init(fan_in, fan_out, constant = 1):
    
    """
 
    目的是合理初始化权重。
 
    参数：
 
    fan_in --行数；
 
    fan_out -- 列数；
 
    constant --常数权重，条件初始化范围的倍数。
 
    return 初始化后的权重tensor.
    """
    low = -constant * np.sqrt(6.0 / (fan_in + fan_out))
 
    high = constant * np.sqrt(6.0 / (fan_in + fan_out))
 
    return tf.random_uniform((fan_in, fan_out),minval = low, maxval = high, dtype = tf.float32)
 
 #%%
 
class AdditiveGaussianNoiseAutoencoder(object):
 

 
    def __init__(self, n_input, n_hidden, transfer_function =tf.nn.sigmoid, optimizer = tf.train.AdamOptimizer(),
     
    scale = 0.1):
     
        self.n_input = n_input
         
        self.n_hidden = n_hidden
         
        self.transfer = transfer_function
         
        self.scale = tf.placeholder(tf.float32)
         
        self.training_scale = scale
         
        network_weights = self._initialize_weights()
         
        self.weights = network_weights
 

 
        self.x = tf.placeholder(tf.float32, [None, self.n_input])
         
        self.hidden = self.transfer(tf.add(tf.matmul(self.x +scale * tf.random_normal((n_input,)),
         
        self.weights['w1']),
         
        self.weights['b1']))
         
        self.reconstruction = tf.add(tf.matmul(self.hidden,self.weights['w2']), self.weights['b2'])
 
 
        self.cost = 0.5 *tf.reduce_sum(tf.pow(tf.subtract(self.reconstruction, self.x), 2.0))
         
        self.optimizer = optimizer.minimize(self.cost)
         
         
        init = tf.global_variables_initializer()
         
        self.sess = tf.Session()
         
        self.sess.run(init)
 
 
    def _initialize_weights(self):
     
        all_weights = dict()
         
        all_weights['w1'] = tf.Variable(xavier_init(self.n_input,self.n_hidden))
         
        all_weights['b1'] = tf.Variable(tf.zeros([self.n_hidden],dtype = tf.float32))
         
        all_weights['w2'] = tf.Variable(tf.zeros([self.n_hidden,self.n_input], dtype = tf.float32))
         
        all_weights['b2'] = tf.Variable(tf.zeros([self.n_input],dtype = tf.float32))
         
        return all_weights
     
     
    def partial_fit(self, X):
     
        cost, opt = self.sess.run((self.cost, self.optimizer),feed_dict = {self.x: X,
         
        self.scale: self.training_scale})
         
        return cost
     
     
    def calc_total_cost(self, X):
     
        return self.sess.run(self.cost, feed_dict = {self.x: X, self.scale:self.training_scale})
     
     
    def transform(self, X):
     
        return self.sess.run(self.hidden, feed_dict = {self.x: X,self.scale:self.training_scale})
     
     
    def generate(self, hidden = None):
     
        if hidden is None:
         
            hidden = np.random.normal(size = self.weights["b1"])
         
        return self.sess.run(self.reconstruction, feed_dict ={self.hidden: hidden})
     
    def reconstruct(self, X):
     
        return self.sess.run(self.reconstruction, feed_dict ={self.x: X,
     
        self.scale: self.training_scale})
     
     
    def getWeights(self): 
     
        return self.sess.run(self.weights['w1'])
     
     
    def getBiases(self): 
     
        return self.sess.run(self.weights['b1'])
 

# In[]。

def get_random_block_from_data(data, batch_size):
 
    start_index = np.random.randint(0, len(data) - batch_size)
     
    return data[start_index:(start_index + batch_size)]


n_samples = len(X)
n_input=len(X[0])
 
#training_epochs = 200
 
batch_size = 128
 
display_step = 1
 
autoencoder = AdditiveGaussianNoiseAutoencoder(n_input,
 
n_hidden,
 
transfer_function =tf.nn.softplus,
 
optimizer =tf.train.AdamOptimizer(learning_rate),
 
scale = 0.01)
 
for epoch in range(training_epochs):
 
    avg_cost = 0.
 
    total_batch = int(n_samples / batch_size)
 
    # Loop over all batches
 
    for i in range(total_batch):
 
        batch_xs = get_random_block_from_data(X, batch_size)
 
 
        # Fit training using batch data
     
        cost = autoencoder.partial_fit(batch_xs)
     
        # Compute average loss
         
        avg_cost += cost / n_samples
     
    # Display logs per epoch step
     
    if epoch % display_step == 0:
     
        print("Epoch:", '%04d' % (epoch + 1), "cost=","{:.9f}".format(avg_cost))
 
# In[] 。

x=[]
x=autoencoder.transform(X)

import pandas as pd

test=pd.DataFrame(index=title,data=x)

output_file=input_file+"_encoder.csv"

test.to_csv(output_file,encoding='gbk') #,header=0