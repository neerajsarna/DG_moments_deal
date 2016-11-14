#!/bin/bash

clear
echo "Testing G120"
BaseDir="xml:/Users/neerajsarna/sciebo/DG_moments_deal/Tests/Test_Reports"

#first we do single core checking
echo "Test Single Core Kn0p1 Odd Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G120_Kn0p1_Odd_SingleCore.xml"
        ./test.out 1 -p ../test_input_files/G120/input_G120_periodic_odd_Kn0p1.in

echo "Test Single Core Kn0p1 Char Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G120_Kn0p1_Char_SingleCore.xml"
        ./test.out 1 -p ../test_input_files/G120/input_G120_periodic_char_Kn0p1.in


#now we do checcking for multicore
# echo "Test Multi Core Kn0p1 Odd Boundary Conditions"
#         export GTEST_OUTPUT="$BaseDir/G120_Kn0p1_Odd_MultiCore.xml"
#         ./test.out 4 -p ../test_input_files/G120/input_G120_periodic_odd_Kn0p1.in

# echo "Test Multi Core Kn0p1 Char Boundary Conditions"
#         export GTEST_OUTPUT="$BaseDir/G120_Kn0p1_Char_MultiCore.xml"
#         ./test.out 4 -p ../test_input_files/G120/input_G120_periodic_char_Kn0p1.in

#Testing for Kn = 0.3
echo "Test Single Core Kn0p3 Odd Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G120_Kn0p3_Odd_SingleCore.xml"
        ./test.out 1 -p ../test_input_files/G120/input_G120_periodic_odd_Kn0p3.in

echo "Test Single Core Kn0p3 Odd Boundary Conditions"
        export GTEST_OUTPUT="$BaseDir/G120_Kn0p3_Char_SingleCore.xml"
        ./test.out 1 -p ../test_input_files/G120/input_G120_periodic_char_Kn0p3.in


#now we do checking for multicore
# echo "Test Multi Core Kn0p3 Odd Boundary Conditions"
#         export GTEST_OUTPUT="$BaseDir/G120_Kn0p3_Odd_MultiCore.xml"
#         ./test.out 4 -p ../test_input_files/G120/input_G120_periodic_odd_Kn0p3.in

# echo "Test Multi Core Kn0p3 Odd Boundary Conditions"
#         export GTEST_OUTPUT="$BaseDir/G120_Kn0p3_Char_MultiCore.xml"
#         ./test.out 4 -p ../test_input_files/G120/input_G120_periodic_char_Kn0p3.in