#!/bin/bash
for file in *.F90
do
sed 's/get_secondary_material(/get_secondary_material(dx,/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
