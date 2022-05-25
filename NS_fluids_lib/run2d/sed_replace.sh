#!/bin/bash
for file in inp*
do
sed 's/ns.continuous_mof=4/ns.continuous_mof=2/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
