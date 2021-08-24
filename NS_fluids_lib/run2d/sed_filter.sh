#!/bin/bash
sed '/[^.0-9 ]/d' $1 > $1_filter_
mv $1_filter_ $1
sed '/^[0-9]/d' $1 > $1_filter_
