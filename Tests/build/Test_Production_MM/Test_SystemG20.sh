#!/bin/bash

clear
echo "Testing G20"
BaseDir=$1

#first we do single core checking
echo "Test Single Core Kn0p1 Odd Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G20_Production_MM_SingleCore.xml"
        ./test.out 1 -p ../test_production_MM/G20/input_G20_Kn0p1.in
