#!/bin/bash

clear
echo "Testing SystemA"
BaseDir="xml:/Users/neerajsarna/sciebo/DG_moments_deal/Tests/Test_Reports"

#first we do single core checking
echo "Test Single Core Kn0p1 Odd Boundary Conditions"
	export GTEST_OUTPUT="$BaseDir/systemA_Kn0p1_Odd_SingelCore.xml"
	./test.out 1 -p ../test_input_files/systemA/input_systemA_ring_odd.in

echo "Test Single Core Kn0p1 Char Boundary Conditions"
	export GTEST_OUTPUT="$BaseDir/systemA_Kn0p1_Char_SingelCore.xml"
	./test.out 1 -p ../test_input_files/systemA/input_systemA_ring_char.in


#now we do checcking for multicore
echo "Test Multi Core Kn0p1 Odd Boundary Conditions"
	export GTEST_OUTPUT="$BaseDir/systemA_Kn0p1_Odd_MultiCore.xml"
	./test.out 4 -p ../test_input_files/systemA/input_systemA_ring_odd.in

echo "Test Multi Core Kn0p1 Char Boundary Conditions"
	export GTEST_OUTPUT="$BaseDir/systemA_Kn0p1_Char_MultiCore.xml"
	./test.out 4 -p ../test_input_files/systemA/input_systemA_ring_char.in 

#Testing for Kn = 1.0
echo "Test Single Core Kn10p0 Odd Boundary Conditions"
	export GTEST_OUTPUT="$BaseDir/systemA_Kn1p0_Odd_SingelCore.xml"
	./test.out 1 -p ../test_input_files/systemA/input_systemA_ring_odd_Kn1p0.in

echo "Test Single Core Kn10p0 Char Boundary Conditions"
	export GTEST_OUTPUT="$BaseDir/systemA_Kn1p0_Char_SingelCore.xml"
	./test.out 1 -p ../test_input_files/systemA/input_systemA_ring_char_Kn1p0.in


#now we do checcking for multicore
echo "Test Multi Core Kn10p0 Odd Boundary Conditions"
	export GTEST_OUTPUT="$BaseDir/systemA_Kn1p0_Odd_MultiCore.xml"
	./test.out 4 -p ../test_input_files/systemA/input_systemA_ring_odd_Kn1p0.in

echo "Test Multi Core Kn10p0 Char Boundary Conditions"
	export GTEST_OUTPUT="$BaseDir/systemA_Kn1p0_Char_MultiCore.xml"
	./test.out 4 -p ../test_input_files/systemA/input_systemA_ring_char_Kn1p0.in 
