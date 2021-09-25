#!/bin/bash
for file in inp*
do
sed 's/constant_viscosity/uncoupled_viscosity/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
