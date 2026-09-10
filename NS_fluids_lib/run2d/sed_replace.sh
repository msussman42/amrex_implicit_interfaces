#!/bin/bash
for file in inp*
do
sed 's/material_extend_velocity/tessellate_elastic_separately/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
