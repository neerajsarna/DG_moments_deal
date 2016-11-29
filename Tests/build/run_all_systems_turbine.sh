#!/bin/sh

echo "Running Tests On All Systems"
clear

#
for VARIABLE in 20 26 35 45 56 71 84 105 120
do 
./run.sh "xml:/home/sarna/DG_moments/Tests/Test_Reports"
done