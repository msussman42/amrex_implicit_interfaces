#!/bin/bash
rm cenx
rm cenx0
rm ceny
rm ceny0
rm linecount
grep "MAT=0 cendir=0 centroid=" run.out > cenx0
grep "MAT=0 cendir=1 centroid=" run.out > ceny0
sed -i 's/TIME=//g' cenx0
sed -i 's/TIME=//g' ceny0
sed -i 's/MAT=0 cendir=0 centroid=//g' cenx0
sed -i 's/MAT=0 cendir=1 centroid=//g' ceny0
grep -c ^ cenx0 > linecount
cat linecount cenx0 > cenx
cat linecount ceny0 > ceny
gfortran convert_to_parametric.f90
./a.out
##set size ratio -1 (gnuplot)

