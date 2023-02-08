#!/bin/bash
sed -i -e -E '/debug_ngrow/s/,[0-9]+);/,local_caller_string);/' Diffusion.cpp
