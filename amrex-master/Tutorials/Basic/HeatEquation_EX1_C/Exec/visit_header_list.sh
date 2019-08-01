#!/bin/bash

#list directory names (Yang Liu via Mehdi Vahab)
ls -w 1 -d plt*/ -tr > Header.visit
#ls -w 1 -d plt*/  > Header.visit
#add header to end of each line
sed -e 's/$/Header/' -i Header.visit
echo new file: Header.visit
