#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
cp ../couplerio.mod .
gfortran solid.f90 ../CouplerIO.f90 -o solid
rm *.bin *.mpf
./solid > solid.out &
../amr3d.Linux.g++.gfortran.ex inputs.whale > run.out 
