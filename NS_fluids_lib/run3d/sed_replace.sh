#!/bin/bash
for file in inp*
do
sed 's/change_max       = 1.1/change_max=1.01/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
