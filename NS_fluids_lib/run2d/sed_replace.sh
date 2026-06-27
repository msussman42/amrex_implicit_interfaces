#!/bin/bash
for file in inp*
do
sed 's/ns\.tiling/ns.ns_tiling/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
