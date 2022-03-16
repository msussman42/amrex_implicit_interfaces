#!/bin/bash
for file in inp*
do
sed 's/ns.visc_abs_tol/mac.visc_abs_tol/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
