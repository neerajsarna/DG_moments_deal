#!/bin/bash

clear
echo "Testing G71"
BaseDir=$1

#first we do single core checking
echo "Test Single Core Kn0p1 Odd Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G71_Production_MM_SingleCore.xml"
        ./test.out 1 -p ../test_production_MM/G71/input_G71_Kn0p1.in
