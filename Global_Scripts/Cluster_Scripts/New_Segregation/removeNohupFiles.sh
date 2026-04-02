#!/bin/sh
# A script to delete the nohup files for all the runs of a particular "architecture" to conserve space. These files are not useful if the simulation ran as expected without errors. All the information is present in the log files too

numberOfMonomers=$1
architecture=$2

numberOfRuns=50

cd b${numberOfMonomers}/${architecture}/
for((i=1; i<=${numberOfRuns}; i++))
do
    rm run$i/nohup.out
done
echo "Removed ${architecture} nohup files!"
