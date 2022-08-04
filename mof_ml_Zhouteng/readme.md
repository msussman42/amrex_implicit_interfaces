# Note

When generating the training data with my subroutine in **generate_data**, 
the domain of the unit cube is [0,1] in each coordinate,
in the legacy code, the domain of the unit is [-0.5,0.5] in each coordinate.

# Compile and run

make compile_gendata: compile the executable for generating training data

make gendata: generate training data and convert to binary

make training: run the example training python file

make compile_gendata: compile the executable for prediction

make gendata: run the example prediction executable

make all: execute all above

make clean: clean the object files and data files

# File list
## generate_data.f90 
Generate training data. 

Calculate the exact centroid from random volume fraction and angles (corresponds with unit normal vector), then get the initial guess from the exact centroid.

At line 11, change the **num_sampling** with a larger number.

All data stored with ascii form.

## convert2binary.py
convert the ascii data to binary form.

## training.py
training with **scikit-learn** library and store the machine learning coefficients.

### inputs and outputs
inputs: initial guess of the angle and volume fractions (3 variables)

outputs: exact angle  (2 variables)

### Parameters for traning

Refer to the scikit-learn manual to set up the training parameters.
In this file, 
neural network algorithm uses default parameters; 
decision tree algorithm set the maximum depth of the tree to 20;
random forest algorithm set the maximum depth of the tree to 20, and count of tree to 10.

[Nueral network](https://scikit-learn.org/stable/modules/generated/sklearn.neural_network.MLPRegressor.html)

[Decision tree](https://scikit-learn.org/stable/modules/generated/sklearn.tree.DecisionTreeRegressor.html#sklearn.tree.DecisionTreeRegressor)

[Random forest](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestRegressor.html)

### Coefficient files for fortran prediction

Neural network: **nn_coef.dat**

Decision tree: **dt_coef.dat**

Random forest: **rf_coef.dat**

The format of the file can refer to **sk2f.py**

## predict.f90

Exanple for prediction.

Read the coefficient file for initialization, and predict with the given input data.

## ml.f90

Modules for Machine learning algorithms