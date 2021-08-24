#!/bin/bash
for file in run*.out
do
grep "isort= 0 im= 2 blob_pressure=" \
  ${file} > ${file}_sumint1_
sed 's/isort= 0 im= 2 blob_pressure=//' \
  ${file}_sumint1_ > ${file}_sumint2_
sed 's/TIME=//' ${file}_sumint2_ > ${file}_pressure_average_
rm ${file}_sumint1_
rm ${file}_sumint2_
done
