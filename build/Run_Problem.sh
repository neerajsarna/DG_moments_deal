#!/bin/bash


# $2 name of the folder where the input files are
# $1 name of the system which has to be tested
for filename in $2/$1/*; do
        ./dg_moments.out 16 -p $filename
done
