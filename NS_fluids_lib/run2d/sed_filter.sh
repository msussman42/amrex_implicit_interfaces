#!/bin/bash
# ^  Matches the starting position within the string
# [^  ] Matches a single character that is not contained within the brackets
sed '/[^.0-9 ]/d' $1 > $1_filter_
mv $1_filter_ $1
sed '/^[0-9]/d' $1 > $1_filter_
mv $1_filter_ $1
