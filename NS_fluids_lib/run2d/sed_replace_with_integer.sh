#!/bin/bash
for file in *.F90
do
sed 's/INTEGER_T/integer/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
