#!/bin/bash
# This script uses the predictConcatenation.c script to find the state of concateation between polymers for all runs of an architecture.

numberOfMonomers=$1
architecture=$2
internalFlag=$3 # The flag indicating whether internal concatenations should be checked for
runIndex=$4 # The index of the run to be launched;
special_simulation=$5 # The name of the initialization procedure
numberOfRuns=50
if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] || [ -z ${internalFlag} ] || [ -z ${runIndex} ] || [ -z ${special_simulation} ]; then
	echo "Not enough arguments passed. Please pass the number of monomer, the architecture name, the internal concatenation flag, the run index, and the initialization procedure as command line arguments."
	exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd $SCRIPT_DIR/../../ && pwd)/"

CHECK_CONCATENATION=${BASE_DIR}Check_Concatenation/
cd ${CHECK_CONCATENATION}Scripts

if [ -z ${special_simulation} ]; then # special simulation not passed
        folderPrefix=b${numberOfMonomers}/
        echo "No special simulation folder passed. Reading files from the default Check_Concatenation/ directory."
else # special simulation passed
        folderPrefix=b${numberOfMonomers}/${special_simulation}/
        echo "Reading files from the special simulation folder '${folderPrefix}'."

fi

rm ../${folderPrefix}${architecture}/concatenations.txt # resetting/removing the earlier concatenations log file
gcc predictConcatenation.c -lm -o con.out -std=gnu99 -I ../../

i=$runIndex
# for((i=1; i<=${numberOfRuns}; i++)) # Can be done for multiple runs if desired
# do
./con.out ${numberOfMonomers} ${architecture} $i ${special_simulation} ${internalFlag} >> ../${folderPrefix}${architecture}/concatenations.txt
# done
