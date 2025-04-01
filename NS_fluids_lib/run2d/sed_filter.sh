#!/bin/bash
# ^  Matches the starting position within the string
# [^  ] Matches a single character that is not contained within the brackets
# using grep if first one or two are digits and third is decimal: 
# grep -E "^[0-9]{1,2}[.]" fname > fname2
# in gnuplot:
# set xrange[xlo:xhi]
# set yrange[ylo:yhi]
sed '/[^.0-9 ]/d' $1 > $1_filter_
mv $1_filter_ $1
sed '/^[0-9]/d' $1 > $1_filter_
mv $1_filter_ $1
