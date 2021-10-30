#!/bin/bash
for file in mat*
do
sed 's/STRANDID= .*/STRANDID=0/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
