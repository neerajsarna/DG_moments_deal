#!/bin/bash

clear
echo "Testing G45"
BaseDir=$1

#first we do single core checking
echo "Test Single Core Kn0p1 Odd Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G45_Production_MM_SingleCore.xml"
        ./test.out 1 -p ../test_production_MM/G45/input_G45_Kn0p1.in
