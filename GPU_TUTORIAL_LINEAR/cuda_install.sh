#!/bin/bash

##first we need to load the anaconda module
module load anaconda/3.8.3

## Now we activate conda
conda init bash

## Now we reload .bashrc
source ~/.bashrc

# Now create the conda environment
conda create -n cuda12env -y
conda activate cuda12env

## Now install CUDA 12.2 Toolkit
conda install -y -c "nvidia/label/cuda-12.2.2" cuda-toolkit

## Now create the folder for local modules
cd ~
mkdir local_modulefiles

## Now append that path to the MODULESPATH variable in .bashrc
echo "export MODULEPATH=$MODULEPATH:${HOME}/local_modulefiles" >> ~/.bashrc

## Now print the cuda-12.2 modulefile to local_modulefiles
echo "#%Module 1.0" > ~/local_modulefiles/cuda-12.2
printf "
#
#  CUDA module for use with 'environment-modules' package:
#
module-whatis  \"Set paths for CUDA-11.1.\"
conflict cuda/10.1
conflict cuda/11.1
set version  12.2
set CUDA_BASE ~/.conda/envs/cuda12env
setenv CUDA_HOME ~/.conda/envs/cuda12env

prepend-path PATH  \${CUDA_BASE}/bin
prepend-path LD_LIBRARY_PATH \${CUDA_BASE}/lib
prepend-path LD_LIBRARY_PATH \${CUDA_BASE}/lib/stubs" >> ${HOME}/local_modulefiles/cuda-12.2

# Lastly, re-run .bashrc again
source ~/.bashrc

