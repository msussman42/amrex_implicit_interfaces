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

exact_f = exact_f.reshape([len(exact_f),1])
inputs = np.hstack((initial_angle,exact_f))
outputs = exact_angle.copy()

nn = MLPRegressor()
nn.fit(inputs, outputs)
nn.output_coef_path = ''

rf = RandomForestRegressor(max_depth=20,n_estimators=10)
rf.fit(inputs, outputs)
rf.output_coef_path = ''

dt = DecisionTreeRegressor(max_depth=20)
dt.fit(inputs, outputs)
dt.output_coef_path = ''

print('R2 score for DT: ',dt.score(inputs[:1000,:],outputs[:1000,:]))
print('R2 score for RF: ',rf.score(inputs[:1000,:],outputs[:1000,:]))
print('R2 score for NN: ',nn.score(inputs[:1000,:],outputs[:1000,:]))

dt_sk2f(dt)
rf_sk2f(rf)
nn_sk2f(nn)