#!/bin/sh

echo "Compiling..........."
mpicxx -I/home/sarna/Documents/Libraries/eigenv3/eigen Test_Eigen.cc -o Test_Eigen
echo "running............."
mpiexec -np 1 ./Test_Eigen
