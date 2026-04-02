#!/bin/bash
# This script runs the Python plotting scripts for the segregation data. It is assumed the data has been downloaded from the cluster.

if [ -z ${POLYMER_BASE} ]; then
	echo "The POLYMER_BASE directory environment variable has not been set!"
	exit 1
fi

NEW_SEGREGATION=${POLYMER_BASE}LAMMPS_runs/cluster-data/new_segregation/

numberOfMonomers=$1
architecture=$2
special_simulation=$3 # optional argument

if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ]
then
	echo "Not enough arguments passed! Please pass the number of monomers and architecture as arguments to the command line."
	echo "An optional argument for the name of a special simulation subfolder can be passed."
	exit 1
fi

cd ${NEW_SEGREGATION}

# Setting up the required directories:
error=0
if ! [ -z ${special_simulation} ]; then
	folderPrefix=b${numberOfMonomers}/Previous_Attempts/${special_simulation}/
	mkdir -p ${folderPrefix}${architecture}/
	echo "Plotting files in the special simulation folder '${folderPrefix}'."
else
	folderPrefix=b${numberOfMonomers}/
	mkdir -p b${numberOfMonomers}/${architecture}/
	echo "No special simulation folder passed. Plotting files in the default new_segregation/b${numberOfMonomers} directory."
fi

mkdir -p ${folderPrefix}${architecture}/Analysis/MonomerDistribution
error=$?
mkdir -p ${folderPrefix}${architecture}/Analysis/TimeSeries
error=$[${error}+$?]
mkdir -p ${folderPrefix}${architecture}/Analysis/Region_CoM_TimeSeries
error=$[${error}+$?]
mkdir -p ${folderPrefix}${architecture}/Analysis/RadialDistribution
error=$[${error}+$?]
mkdir -p ${folderPrefix}${architecture}/Analysis/MinDistance
error=$[${error}+$?]

if ! [ ${error} -eq 0 ]; then
	echo "Something went wrong while making the required directories. Please try again."
	exit 1
fi

cd ${NEW_SEGREGATION}Analysis/

# Performing analysis: Plotting monomer density 
python plotMonomerDensity.py ${numberOfMonomers} ${architecture}

# Plotting time series/CoM distance:
python plotCoMDistribution.py ${numberOfMonomers} ${architecture}
error=$?
if ! [ ${error} -eq 0 ]; then
	echo "Something went wrong while plotting the CoM Time Series. Terminating."
	exit 1
fi

# Plotting minimum distances (infinite cylinder):
if [[ ${special_simulation} == inf_* ]]; then
	echo "Plotting Minimum Distances for ${special_simulation}"
	python plotMinimumDistance.py ${numberOfMonomers} ${architecture}
	error=$?
	if ! [ ${error} -eq 0 ]; then
		echo "Something went wrong while plotting the minimum distances. Terminating."
		exit 1
	fi
else
	echo "Not Plotting Minimum Distances for ${special_simulation}."
fi

# Plotting distribution of segregation time:
python plotSegTimeDistribution.py ${numberOfMonomers} ${architecture} both
error=$?

if ! [ ${error} -eq 0 ]; then
	echo "Something went wrong while plotting the segregation time distribution."
fi


echo "Done!"
