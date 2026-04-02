#!/bin/sh
# This script should be run after the CoM correction has been done, if required. This script deletes the original CoM files and replaces them by the corrected ones (by renaming the corrected files to the original files)
architecture=$1 # setting the first argument passed as the architecture
numberOfRuns=50
numberOfMonomers=$2
echo "Architecture = ${architecture}"
cd b${numberOfMonomers}/${architecture}/
for((i=1; i<=${numberOfRuns}; i++))
#i=1
do
    cd run$i/
    rm com*
    for((j=1; j<=2; j++))
    do
	mv corrected_com$j.dat com$j.dat
    done
    cd ..
done
