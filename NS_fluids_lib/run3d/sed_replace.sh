#!/bin/bash
for file in inp*
do
sed 's/file_name_digits=7/file_name_digits=8/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
