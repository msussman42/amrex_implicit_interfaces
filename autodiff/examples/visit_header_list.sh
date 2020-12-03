#!/bin/bash

#list directory names (Yang Liu via Mehdi Vahab)
ls -1 plt*/Header | tee Header.visit
#add header to end of each line
echo new file: Header.visit
