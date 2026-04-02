#!/bin/bash
# This is a script that runs execute_calculation.sh using nohup

# input parameters:
numberOfMonomers=$1
architecture=$2


if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] # -z flag checks whether the variable has zero length
then
        echo "Not enough arguments passed! Please pass number of monomers and the architecture while executing this script!"
        exit 1
fi

# Making folder
mkdir -p b${numberOfMonomers}/${architecture}/
nohup bash execute_calculation.sh ${numberOfMonomers} ${architecture} > b${numberOfMonomers}/${architecture}/nohup.out &
sleep 1
