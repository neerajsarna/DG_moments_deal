#!/bin/bash

for filename in ../$2/$1/*; do
        ./dg_moments.out 16 -p $filename
done