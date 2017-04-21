#!/bin/sh

for i in {5..16}
do
	./test.out 1 -p ../../Input_Files/1D_Systems/input_N$i.in -p ../../Input_Files/Input_Testing/Input_Inflow_1D/input_poisson_1D_meshworker.in 
done
