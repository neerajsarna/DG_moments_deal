#!/bin/bash


# how to run the skript
# 1st argument is the system to be tested (with the G)
# 2nd argument is the name of the test case which has to be studied
clear
echo "CMAKE"

export OMP_NUM_THREADS=16
rm dg_moments.out
rm -rf CMakeFiles/ CMakeCache.txt
cmake ../.
make release


./Run_Problem.sh $1 $2
# read -p "Continue with execution?" yn
#     case $yn in
#         [Yy]* ) ./Test_SystemG$1.sh $2; break;;
#         [Nn]* ) exit;;
#         * ) echo "Please answer yes or no.";;
#     esac
