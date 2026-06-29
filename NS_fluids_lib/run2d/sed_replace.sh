#!/bin/bash
for file in inp*
do
sed 's/LSA_nsteps_power_method/LSA_nsteps_krylov_subspace_method/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
