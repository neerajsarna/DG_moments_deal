#!/bin/bash

# argument-1: is the type of boundary condition
# argument-2: the value of the Knudsen number


for filename in ../Input_Poisson_Heat_BGK/GA/$2/$1/*; do
        ./dg_moments.out 16 -p $filename
done
