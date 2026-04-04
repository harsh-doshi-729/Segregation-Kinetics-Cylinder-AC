#!/bin/bash
# This script uses the predictConcatenation.c script to find the state of concateation between polymers for all runs of an architecture.

numberOfMonomers=$1
architecture=$2
internalFlag=$3 # The flag indicating whether internal concatenations should be checked for
special_simulation=$4 # optional argument
numberOfRuns=50
if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] || [ -z ${internalFlag} ]; then
	echo "Not enough arguments passed. Please pass the number of monomer, the architecture name, and the internal concatenation flag as command line arguments."
	echo "An optional argument for the special simulation identifier can be passed as well to store the files in a special directory"
	exit 1
fi

BASE_DIR=$(dirname $(dirname $(dirname $(readlink -f $0))))
echo "Base directory for the Check_Concatenation scripts: ${BASE_DIR}"
CHECK_CONCATENATION=${BASE_DIR}/Check_Concatenation/
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

for((i=1; i<=${numberOfRuns}; i++))
do
	./con.out ${numberOfMonomers} ${architecture} $i ${internalFlag} >> ../${folderPrefix}${architecture}/concatenations.txt
done
