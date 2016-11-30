#!/bin/bash

clear
echo "Testing G56"
BaseDir=$1

#Testing for Kn = 0.1
echo "Test Single Core Kn0p1 Odd Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G56_Kn0p1_Odd_SingleCore.xml"
        ./test.out 1 -p ../test_poisson_heat_bgk/G56/input_G56_periodic_odd_Kn0p1.in

echo "Test Single Core Kn0p1 Char Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G56_Kn0p1_Char_SingleCore.xml"
        ./test.out 1 -p ../test_poisson_heat_bgk/G56/input_G56_periodic_char_Kn0p1.in



#Testing for Kn = 0.3
echo "Test Single Core Kn0p3 Odd Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G56_Kn0p3_Odd_SingleCore.xml"
        ./test.out 1 -p ../test_poisson_heat_bgk/G56/input_G56_periodic_odd_Kn0p3.in

echo "Test Single Core Kn0p3 Odd Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G56_Kn0p3_Char_SingleCore.xml"
        ./test.out 1 -p ../test_poisson_heat_bgk/G56/input_G56_periodic_char_Kn0p3.in

