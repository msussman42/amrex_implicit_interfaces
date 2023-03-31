#!/bin/bash
for file in inp*
do
sed 's/2 = CMOF/1=CMOF/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
