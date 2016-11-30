#!/bin/bash

clear
echo "Testing G105"
BaseDir=$1

#first we do single core checking
echo "Test Single Core Kn0p1 Odd Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G105_Production_MM_SingleCore.xml"
        ./test.out 1 -p ../test_production_MM/G105/input_G105_Kn0p1.in
