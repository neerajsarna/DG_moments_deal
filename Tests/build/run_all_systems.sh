#!/bin/sh

# A file for performing tests on all the systems that have been loaded
# First argument is the test which has to be run

echo "Running Tests On All Systems"
clear

#
for VARIABLE in A 20 26 35 45 56 71 84 105 120
do 
./run.sh $VARIABLE $1
done
