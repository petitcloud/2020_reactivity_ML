import sys, os, time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
import numpy as np
from math import log10
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

print('\n\nML Nucleophilicity prediction code, Made by B.LEE')
print('Trained by 726 molecule. RMSE = 3.26')
print(''' 
Your file should be in this format,
------------------------------------------------------
4-aminopyridine  NC1=CC=NC=C1  3.04959102  0.45322  0.30205  0.264    0.26396  0.00329  0.00329  0.00054 -0.00799 -0.00799 -0.01196 -0.01196 -0.12522 -0.12522 
glycineamide   NCC(N)=O   2.847399712  0.89849  0.064    0.06399  0.03633 -0.00104 -0.00112 -0.00352 -0.00428 -0.00825 -0.02198 -0.02261
...
------------------------------------------------------
  (Name)      (SMILES)   (Adiabatic I.E.)  (Parr function, sorted in descending order)

Name -- Name of the molecule
SMILES -- SMILES expression of molecular structure
Adiabatic I.E. -- Ionization energy calculated from fully relaxed (n-1) electron structure
Parr function -- Spin density of fully relaxed (n-1) electron structure, sorted in descending order. 
                 ~30 Parr function values can be provided. write down to the 30th value if the number of atom in molecule is larger than 30.


''')
print('Type the path of list file !')

dirc = input()

#dirc = sys.argv[1]

# model parameters
N_LAYERS = 2
WEIGHTS = [400, 400]
INPUT_SIZE = 30
BATCH_SIZE = 100000
LEARNING_RATE = 1e-03
LAMBDA = 1e-03
DROPOUT_RATE = 0.8
EPOCH = 100001


with open(dirc, 'r') as inputfile:
 indexline = inputfile.readlines()

mol_name = []
mol_smiles = []
mol_indices = []
mol_parr = np.zeros((len(indexline),INPUT_SIZE))

for i in range(len(indexline)):
 mol_name.append(indexline[i].split()[0])
 mol_smiles.append(indexline[i].split()[1])
 mol_indices.append(float(indexline[i].split()[2]))
 temp_list = indexline[i].split()[3:]
 for k in range(min(INPUT_SIZE, len(temp_list))):
  mol_parr[i][k] = pow(10, mol_indices[-1]*float(temp_list[k]))

fp_mat = []
for line in mol_smiles:
 rdtemp = Chem.MolFromSmiles(line)
 rd_mor = MACCSkeys.GenMACCSKeys(rdtemp)
 fp_mat.append(np.array([int(x) for x in rd_mor.ToBitString()]))

if not N_LAYERS == len(WEIGHTS):
 print('The number of layers and Dimension of weights are not same!')
 print('The number of layers: {}, Dimension of weights: {}' .format(N_LAYERS, len(WEIGHTS)))
 sys.exit()

# Import data using gen_fingerprint_array and gen_parr_array function

if not len(fp_mat) == len(mol_parr):
 print('The number of features are not same !!')
 print('The number of fingerprint features: {}, the number of parr function features: {}' .format(fp_mat.shape[0], parr_ind_mat.shape[0]))
 sys.exit()

fp_mat = np.array(fp_mat)
mol_parr = np.array(mol_parr)

# Integrating two features to single array
feature = np.zeros((fp_mat.shape[0], fp_mat.shape[1]+INPUT_SIZE))
feature[:,:fp_mat.shape[1]] = fp_mat
feature[:,fp_mat.shape[1]:] = mol_parr[:,:INPUT_SIZE]
feature = tf.constant(feature, tf.float64)

# Build dataset
train_dataset = tf.data.Dataset.from_tensor_slices(feature).batch(BATCH_SIZE)
iterator_dataset = train_dataset.make_initializable_iterator()
iterator_dataset.initializer
next_element = iterator_dataset.get_next()

# Linear model with n hidden layers
WEIGHTS = [int(next_element.shape[1])] + WEIGHTS
W_hidden = []
b_hidden = []
for k in range(N_LAYERS):
 W_hidden.append(tf.Variable(tf.random_uniform([WEIGHTS[k],WEIGHTS[k+1]], 0.0, 0.01, dtype=tf.float64), name='W_h{}' .format(k+1)))
 b_hidden.append(tf.Variable(tf.zeros([WEIGHTS[k+1]], dtype=tf.float64), name='b_h{}' .format(k+1)))

W_o = tf.Variable(tf.random_uniform([WEIGHTS[-1],1], 0.0, 0.01, dtype=tf.float64), name='W_o')
b_o = tf.Variable(tf.zeros([1], dtype=tf.float64), name='b_o')

saver = tf.train.Saver([W_hidden[i] for i in range(N_LAYERS)] + [b_hidden[i] for i in range(N_LAYERS)]  +[W_o, b_o])

def model(x, DROPOUT_RATE):
 hidden_layer = []
 for l in range(N_LAYERS):
  if l == 0:
   hidden_layer.append(tf.nn.dropout(tf.nn.relu(tf.matmul(x, W_hidden[l])+b_hidden[l], name='hidden{}' .format(l+1)), DROPOUT_RATE))
  else:
   hidden_layer.append(tf.nn.dropout(tf.nn.relu(tf.matmul(hidden_layer[l-1], W_hidden[l])+b_hidden[l], name='hidden{}' .format(l+1)), DROPOUT_RATE))
 y_= tf.transpose(tf.matmul(hidden_layer[-1], W_o) + b_o)
 return y_

# Define loss and trainset
def loss(x, y, DROPOUT_RATE):
 y_ = model(x, DROPOUT_RATE)
 return tf.reduce_mean(tf.square(y - y_))

def rsquare(y, lo):
 ss = sess.run(tf.reduce_sum(tf.square(y-tf.reduce_mean(y))))
 return (1-lo*int(y.shape[0])/ss)

now1 = time.time()
# Prediction step
with tf.Session() as sess:
 sess.run(iterator_dataset.initializer)
 saver.restore(sess, './train_ckpt/trained_weights.ckpt')
 x = sess.run(next_element) # Instead of using gen_sim_data
 y = sess.run(model(x,1))

now2 = time.time()
print("Calculation time: {} seconds" .format(now2-now1))

if __name__ == '__main__':
 print('----------------------')
 print('<folder_name>\t<prediced_nucleophilicity>')
 for i in range(len(mol_name)):
  print(mol_name[i], float(y[:,i]))
