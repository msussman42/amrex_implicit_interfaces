#!/bin/bash
num=0
for filedir in ./data*
do
cd ${filedir}
for file in run*.out
do
grep "isort= 0 im= 2 blob_pressure=" \
  ${file} > ${file}_sumint1_
sed 's/isort= 0 im= 2 blob_pressure=//' \
  ${file}_sumint1_ > ${file}_sumint2_
newfile=_sumint_${file}_sumint_
sed 's/TIME=//' ${file}_sumint2_ > ${newfile}
rm ${file}_sumint1_
rm ${file}_sumint2_
newfile_up_one="$newfile$num"
echo ${newfile}
echo ${newfile_up_one}
mv ${newfile} ../${newfile_up_one}
num=$((num+1))
done
cd ..
done
ls ./_sumint_*
cat _sumint_* >> pressure_average_data
