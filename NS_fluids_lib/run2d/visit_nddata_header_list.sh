#!/bin/bash

ls -w -1 -tr ndd* | grep -o "nddataPLT[0-9]*" > header0.visit
echo temp file: header0.visit
sed 's/$/\/Header/' header0.visit > header.visit
rm header0.visit
echo file: header.visit

