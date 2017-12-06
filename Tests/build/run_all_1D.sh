#!/bin/sh

for i in {33..50}
do
        rem=$(($i % 2))
        if [ "$rem" -ne "0" ]; then
		./test.out 1 -p ../../Input_Files/1D_Systems/input_N$i.in -p ../../Input_Files/Input_Inflow_1D/input_poisson_1D_meshworker.in
                echo "<<<<<<<<<<<<<<<Solving For $i>>>>>>>>>>>>>>>>>"
        fi
done
