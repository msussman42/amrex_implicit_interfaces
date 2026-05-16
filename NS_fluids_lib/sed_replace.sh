#!/bin/bash
for file in *.F90
do
sed 's/VOFTOL/VOFTOL_MATERIAL/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
