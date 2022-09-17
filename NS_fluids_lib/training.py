from numpy.random import seed
seed(5)
import numpy as np
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import explained_variance_score
from sk2f import *

data_dir = '.'

exact_f = np.load(data_dir+'/exact_f.npy')
exact_angle = np.load(data_dir+'/exact_angle.npy')
exact_centroid = np.load(data_dir+'/exact_centroid.npy')
initial_angle = np.load(data_dir+'/initial_angle.npy')

exact_f_shape=exact_f.shape
print ('exact_f_shape ',exact_f_shape)
exact_angle_shape=exact_angle.shape
print ('exact_angle_shape ',exact_angle_shape)
initial_angle_shape=initial_angle.shape
print ('initial_angle_shape ',initial_angle_shape)

num_sampling=exact_angle_shape[0]
print ('num_sampling ',num_sampling)
num_shape_dim=len(exact_angle_shape)

amrex_spacedim=0
if (num_shape_dim==1):
    amrex_spacedim=2
if (num_shape_dim==2):
    amrex_spacedim=exact_angle_shape[1]+1
print ('amrex_spacedim (python) ',amrex_spacedim)

exact_f = exact_f.reshape([num_sampling,1])
initial_angle=initial_angle.reshape([num_sampling,amrex_spacedim-1])
exact_angle=exact_angle.reshape([num_sampling,amrex_spacedim-1])
inputs = np.hstack((initial_angle,exact_f))

outputs = exact_angle.copy()
if (amrex_spacedim==2):
    outputs=outputs.reshape([num_sampling, ])
if (amrex_spacedim==3):
    outputs=outputs.reshape([num_sampling,amrex_spacedim-1])

inputs_shape=inputs.shape
outputs_shape=outputs.shape

print ('inputs_shape ',inputs_shape)
print ('outputs_shape ',outputs_shape)

nn = MLPRegressor()
nn.fit(inputs, outputs)
nn.output_coef_path = ''

rf = RandomForestRegressor(max_depth=20,n_estimators=10)
rf.fit(inputs, outputs)
rf.output_coef_path = ''

#dt = DecisionTreeRegressor(max_depth=20)
dt = DecisionTreeRegressor(max_depth=40)
dt.fit(inputs, outputs)
dt.output_coef_path = ''

outputs=outputs.reshape([num_sampling,amrex_spacedim-1])

print('R2 score for DT: ',dt.score(inputs[:1000,:],outputs[:1000,:]))
print('R2 score for RF: ',rf.score(inputs[:1000,:],outputs[:1000,:]))
print('R2 score for NN: ',nn.score(inputs[:1000,:],outputs[:1000,:]))


print('\n')
print('For input vector (init_angle1,init_angle2,f)',inputs[0,:], '\n predictions are:')
print('NN: ',nn.predict(inputs[0,:].reshape(1,-1)))
print('DT: ',dt.predict(inputs[0,:].reshape(1,-1)))
print('RF: ',rf.predict(inputs[0,:].reshape(1,-1)))
print('Expected output: ',outputs[0,:].reshape(1,-1))

dt_sk2f(dt)
rf_sk2f(rf)
nn_sk2f(nn)
