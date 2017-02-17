#!/bin/bash


clear
echo "CMAKE"

export OMP_NUM_THREADS=16
rm dg_moments.out
rm -rf CMakeFiles/ CMakeCache.txt
cmake ../.
make release

#./run_systemA.sh char Kn0p1
#./run_systemA.sh odd Kn0p1

#./run_systemA.sh char Kn1p0
./run_systemA.sh odd Kn1p0

#./run_systemA.sh char Kn10p0
#./run_systemA.sh odd Kn10p0



