#!/bin/bash


# how to run the skript
# 1st argument is the system to be tested
# 2nd argument is the name of the test case which has to be studied
clear
echo "CMAKE"

BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "Base Directory of the present script "
echo $BASEDIR

rm -rf CMakeFiles/ CMakeCache.txt
cmake -D G$1=ON ../.
make debug


./$2/Test_System.sh "xml:$BASEDIR/$2/Test_ReportsG$1" $1
# read -p "Continue with execution?" yn
#     case $yn in
#         [Yy]* ) ./Test_SystemG$1.sh $2; break;;
#         [Nn]* ) exit;;
#         * ) echo "Please answer yes or no.";;
#     esac
