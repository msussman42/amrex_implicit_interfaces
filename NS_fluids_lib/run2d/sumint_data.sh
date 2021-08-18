#!/bin/bash
for file in run*.out
do
grep "user_comp2 (1..ncomp_sum_int_user2) 1 sum_int_user2" \
  ${file} > ${file}_sumint1_
sed 's/user_comp2 (1..ncomp_sum_int_user2) 1 sum_int_user2//' \
  ${file}_sumint1_ > ${file}_sumint2_
sed 's/TIME=//' ${file}_sumint2_ > ${file}_sumint_
rm ${file}_sumint1_
rm ${file}_sumint2_
done
