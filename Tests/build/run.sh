#!/bin/bash

clear
echo "Cmake"
rm -rf CMakeFiles/ CMakeCache.txt
cmake -D G$1=ON ../.
make debug


./Test_SystemG$1.sh $2
# read -p "Continue with execution?" yn
#     case $yn in
#         [Yy]* ) ./Test_SystemG$1.sh $2; break;;
#         [Nn]* ) exit;;
#         * ) echo "Please answer yes or no.";;
#     esac
