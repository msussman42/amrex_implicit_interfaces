#!/bin/bash

ls -w -1 -tr ndd* | grep -o "nddataPLT[0-9]*" APPEND /Header > header.visit
echo new file: header.visit
