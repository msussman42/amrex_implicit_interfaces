#!/bin/bash
for file in *.F90
do
sed 's/irz.eq.3/irz.eq.COORDSYS_CYLINDRICAL/' ${file} > ${file}_mod_
cp ${file}_mod_ ${file}
rm ${file}_mod_
done
