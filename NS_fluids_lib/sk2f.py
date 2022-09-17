import numpy as np
def sk2f(write_coef):
    if (write_coef.type == 'Neural_Network'):
        nn_sk2f(write_coef)
    elif (write_coef.type == 'Random_Forest'):
        rf_sk2f(write_coef)
    elif (write_coef.type == 'Decision_Tree'):
        dt_sk2f(write_coef)
    else:
        print('wrong type')

def nn_sk2f(nn):

    file_name = 'nn_coef.dat'
    nn_output  = open(nn.output_coef_path+file_name, 'w')
    nn_output.write('! n_layers\n%d\n'%(nn.n_layers_))
    nn_output.write('! layer_sizes\n')
    nn_output.write('%d \n'%np.shape(nn.coefs_[0])[0])
    for i in range(len(nn.hidden_layer_sizes)):
        nn_output.write('%d '%nn.hidden_layer_sizes[i])
    nn_output.write('%d\n'%nn.n_outputs_)
    nn_output.write('! intercepts\n')

    for i in range(len(nn.intercepts_)):
        for j in range(len(nn.intercepts_[i])):
            nn_output.write('%f '%nn.intercepts_[i][j])
        nn_output.write('\n')
    nn_output.write('! coefs\n')
    for i in range(len(nn.coefs_)):
        nn_output.write('!! layer%d\n'%i)
        coef = np.transpose(nn.coefs_[i])
        for j in range(len(coef)):
            for k in range(len(coef[j])):
                nn_output.write('%f '%coef[j][k])
            nn_output.write('\n')
    nn_output.write('! activations\n')
    nn_output.write(nn.activation)
    nn_output.write('\n')
    nn_output.write('! out_activations\n')
    nn_output.write(nn.out_activation_)
    nn_output.write('\n')

def dt_sk2f(dt):

    dt_output = open(dt.output_coef_path+'dt_coef.dat','w')

    dt_output.write('! node_count\n')
    dt_output.write('%d\n'%dt.tree_.node_count)
    dt_output.write('! n_features\n')
    dt_output.write('%d\n'%dt.tree_.n_features)
    dt_output.write('! n_outputs\n')
    dt_output.write('%d\n'%dt.tree_.n_outputs)
    dt_output.write('! max_depth\n')
    dt_output.write('%d\n'%dt.tree_.max_depth)

    for i in range(dt.tree_.node_count):
        dt_output.write('! node %d\n'%(i+1))
        dt_output.write('%d\n'%(dt.tree_.children_left[i]+1))
        dt_output.write('%d\n'%(dt.tree_.children_right[i]+1))
        dt_output.write('%d\n'%(dt.tree_.feature[i]+1))
        dt_output.write('%f\n'%dt.tree_.threshold[i])
        for j in range(dt.tree_.n_outputs):
            dt_output.write('%f '%dt.tree_.value[i,j])
        dt_output.write('\n')

def rf_sk2f(rf):

    rf_output = open(rf.output_coef_path+'rf_coef.dat','w')

    rf_output.write('! tree_count\n')
    rf_output.write('%d\n'%len(rf.estimators_))

    trees = rf.estimators_

    for j in range(len(rf.estimators_)):
        tree1 = trees[j].tree_
        rf_output.write('! node_count\n')
        rf_output.write('%d\n'%tree1.node_count)
        rf_output.write('! n_features\n')
        rf_output.write('%d\n'%tree1.n_features)
        rf_output.write('! n_outputs\n')
        rf_output.write('%d\n'%tree1.n_outputs)
        rf_output.write('! max_depth\n')
        rf_output.write('%d\n'%tree1.max_depth)

        for i in range(tree1.node_count):
            rf_output.write('! node %d\n'%(i+1))
            rf_output.write('%d\n'%(tree1.children_left[i]+1))
            rf_output.write('%d\n'%(tree1.children_right[i]+1))
            rf_output.write('%d\n'%(tree1.feature[i]+1))
            rf_output.write('%f\n'%tree1.threshold[i])
            for j in range(tree1.n_outputs):
                rf_output.write('%f '%tree1.value[i,j])
            rf_output.write('\n')
