#!/bin/bash

clear
echo "Testing $2"
BaseDir=$1
export count=0
echo "Value of count is"
echo $count
for filename in ../test_heated_cavity_MM/G$2/*; do
         export filetowrite="InputFile${count}"
         export GTEST_OUTPUT="$BaseDir/$filetowrite.xml"
         echo "Writting Test Reports to\n"
         echo $GTEST_OUTPUT
         count=$((count+1))
        ./test.out 16 -p $filename
done
