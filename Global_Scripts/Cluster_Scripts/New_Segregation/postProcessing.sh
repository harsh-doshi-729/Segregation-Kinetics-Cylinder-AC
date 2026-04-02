#!/bin/bash
# This bash script contains the commands to be run in order to perform the post processing for the segregated runs.
# This script is invoked in execute_postProcessing.sh

numberOfMonomers=$1
architecture=$2
run=$3
axisLength=$4
destinationLabel=$5



# Assuming that the C files have already been compiled

# Calculating monomer density (region wise):
echo "Calculating region wise monomer distribution:"
./ge.out ${numberOfMonomers} ${architecture} ${run} ${destinationLabel} 0 ${axisLength}; 
error=$?
if ! [ ${error} -eq 0 ]; then
	echo "Something went wrong while calculating the monomer density."
fi
echo ""

# Calculating radial monomer density:
echo "Calculating radial monomer density:"
./radial.out ${numberOfMonomers} ${architecture} ${run} 0.1 # binwidth is taken as 0.1
error=$?
if ! [ ${error} -eq 0 ]; then
	echo "Something went wrong while calculating the radial monomer density."
fi
echo ""

# Calculating the CoM time series region wise:
echo "Calculating region wise CoM Time Series:"
./com.out ${numberOfMonomers} ${architecture} ${run} ${destinationLabel}
error=$?
if ! [ ${error} -eq 0 ]; then
	echo "Something went wrong while calculating the regional CoM time series."
fi
