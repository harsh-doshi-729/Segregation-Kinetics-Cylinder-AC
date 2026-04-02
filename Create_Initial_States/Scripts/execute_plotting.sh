#!/bin/bash
# This script plots the various quantities that have been calculated from the simulation

numberOfMonomers=$1
architecture=$2
# special_simulation=$3 # Optional argument

if [ -z ${POLYMER_BASE} ]; then
	echo "The POLYMER_BASE directory environment variable has not been set!"
	exit 1
fi
CREATE_INITIAL_STATES=${POLYMER_BASE}LAMMPS_runs/cluster-data/new_segregation/Create_Initial_States/
cd ${CREATE_INITIAL_STATES}

if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ]; then
	echo "Please pass the number of monomers and architecture as command line arguments!"
	# echo "An optional argument for the special simulation identifier can be passed as well to store the files in a special directory"
	exit 1
fi

# # Going to the special simulation directory if passed:
# if ! [ -z ${special_simulation} ]; then
# 	folderPrefix=b${numberOfMonomers}/Previous_Attempts/${special_simulation}/
# 	echo "Plotting files in the special simulation folder '${folderPrefix}'."
# else
# 	folderPrefix=b${numberOfMonomers}/
# 	echo "No special simulation folder passed. Plotting files in the default Create_Initial_States/${folderPrefix} directory."
# fi

# Plotting the distributions:
cd Scripts/
python plotMonomerDensity.py ${numberOfMonomers} ${architecture}
error=$?
if ! [ ${error} -eq 0 ]; then
        echo "Something went wrong while plotting the monomer density."
fi

# python plotBondEnergyDistribution.py ${numberOfMonomers} ${architecture}

# error=$?
# if ! [ ${error} -eq 0 ]
# then
#         echo "Something went wrong while plotting the bond energy distribution."
# fi

# Plotting single snapshot distributions:
binWidth=1 # bin width for the single snapshot distribution
python plotSingleSnapshotDistribution.py ${numberOfMonomers} ${architecture} ${binWidth} # For all runs
error=$?
if ! [ ${error} -eq 0 ]
then
        echo "Something went wrong while plotting the individual single-snapshot monomer distributions."
fi
python plotSingleSnapshotDistribution.py ${numberOfMonomers} ${architecture} ${binWidth} 0 # Plotting average
error=$?
if ! [ ${error} -eq 0 ]
then
        echo "Something went wrong while plotting the average single-snapshot monomer distribution."
fi

# Plotting the time series of the polymer and regions:
python plotCoMTimeSeries.py ${numberOfMonomers} ${architecture}
error=$?
if ! [ ${error} -eq 0 ]
then
        echo "Something went wrong while plotting the CoM time series for each run."
fi

# Plotting pair-correlation functions:
binWidth=0.1 # for the pair_correlation
python plotPairCorrelation.py ${numberOfMonomers} ${architecture} ${binWidth} monomer 1 0
error=$?
if ! [ ${error} -eq 0 ]
then
        echo "Something went wrong while plotting the pair correlation for each run."
fi
