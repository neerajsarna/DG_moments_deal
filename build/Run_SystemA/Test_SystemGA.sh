#!/bin/bash

clear
echo "Running SystemA for A Ring"

echo "Char Kn = 0.1"
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_char_Kn0p1_degree_1.in 
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_char_Kn0p1_degree_2.in
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_char_Kn0p1_degree_3.in
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_char_Kn0p1_degree_4.in

echo "Char Kn = 1.0"
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_char_Kn1p0_degree_1.in 
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_char_Kn1p0_degree_2.in
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_char_Kn1p0_degree_3.in
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_char_Kn1p0_degree_4.in


echo "Char Kn = 10.0"
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_char_Kn10p0_degree_1.in 
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_char_Kn10p0_degree_2.in
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_char_Kn10p0_degree_3.in
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_char_Kn10p0_degree_4.in


echo "Odd Kn = 0.1"
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_odd_Kn0p1_degree_1.in 
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_odd_Kn0p1_degree_2.in
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_odd_Kn0p1_degree_3.in
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_odd_Kn0p1_degree_4.in

echo "Odd Kn = 1.0"
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_odd_Kn1p0_degree_1.in 
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_odd_Kn1p0_degree_2.in
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_odd_Kn1p0_degree_3.in
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_odd_Kn1p0_degree_4.in


echo "Odd Kn = 10.0"
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_odd_Kn10p0_degree_1.in 
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_odd_Kn10p0_degree_2.in
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_odd_Kn10p0_degree_3.in
	./dg_moments.out 16 -p ../Input_Poisson_Heat_BGK/systemA/input_systemA_ring_odd_Kn10p0_degree_4.in
