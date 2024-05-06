#!/bin/bash
for file in inp*
do
sed 's/MOFITERMAX=15/MOFITERMAX=30/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
