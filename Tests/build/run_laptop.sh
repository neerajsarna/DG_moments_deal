#!/bin/bash

clear
echo "Cmake"
rm -rf CMakeFiles/ CMakeCache.txt
cmake ../

make release


read -p "Continue with execution?" yn
    case $yn in
        [Yy]* ) ./Test_SystemG$1.sh; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
