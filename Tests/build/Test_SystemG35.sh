#!/bin/bash

clear
echo "Testing G35"
BaseDir=$1

#first we do single core checking
echo "Test Single Core Kn0p1 Odd Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G35_Kn0p1_Odd_SingleCore.xml"
        ./test.out 1 -p ../test_input_files/G35/input_G35_periodic_odd_Kn0p1.in

echo "Test Single Core Kn0p1 Char Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G35_Kn0p1_Char_SingleCore.xml"
        ./test.out 1 -p ../test_input_files/G35/input_G35_periodic_char_Kn0p1.in



#Testing for Kn = 0.3
echo "Test Single Core Kn0p3 Odd Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G35_Kn0p3_Odd_SingleCore.xml"
        ./test.out 1 -p ../test_input_files/G35/input_G35_periodic_odd_Kn0p3.in

echo "Test Single Core Kn0p3 Odd Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G35_Kn0p3_Char_SingleCore.xml"
        ./test.out 1 -p ../test_input_files/G35/input_G35_periodic_char_Kn0p3.in

