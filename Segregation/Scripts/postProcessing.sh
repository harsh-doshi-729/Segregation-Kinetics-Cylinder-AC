#!/bin/bash
# This bash script contains the commands to be run in order to perform the post processing for the segregated runs.
# This script is invoked in execute_postProcessing.sh

numberOfMonomers=$1
architecture=$2
run=$3
special_simulation=$4 # The name of the initialization procedure
destinationLabel=$5
axisLength=$6 # Optional argument for finite cylinder simulations

if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] || [ -z ${run} ] || [ -z ${special_simulation} ] || [ -z ${destinationLabel} ]; then
	echo "Not enough arguments passed! Please pass the number of monomers, the architecture, the run index, the name of the initialization procedure, and the destination label as arguments while invoking the script"
	echo "For example: $ postProcessing.sh 200 Arc-0 1 fene_recenter Segregation"
	echo "An optional argument for the axis length can also be passed for finite cylinder simulations. For example: $ postProcessing.sh 200 Arc-0 1 fene_recenter Segregation 31.2"
	exit 1
fi

# Assuming that the C files have already been compiled

# Calculating monomer density (region wise):
echo "Calculating region wise monomer distribution:"
./ge.out ${numberOfMonomers} ${architecture} ${run} ${destinationLabel} ${special_simulation} 0 ${axisLength}; 
error=$?
if ! [ ${error} -eq 0 ]; then
	echo "Something went wrong while calculating the monomer density."
fi
echo ""

# Calculating radial monomer density:
echo "Calculating radial monomer density:"
./radial.out ${numberOfMonomers} ${architecture} ${run} 0.1 ${special_simulation} # binwidth is taken as 0.1
error=$?
if ! [ ${error} -eq 0 ]; then
	echo "Something went wrong while calculating the radial monomer density."
fi
echo ""

# Calculating the CoM time series region wise:
echo "Calculating region wise CoM Time Series:"
./com.out ${numberOfMonomers} ${architecture} ${run} ${special_simulation} ${destinationLabel}
error=$?
if ! [ ${error} -eq 0 ]; then
	echo "Something went wrong while calculating the regional CoM time series."
fi
