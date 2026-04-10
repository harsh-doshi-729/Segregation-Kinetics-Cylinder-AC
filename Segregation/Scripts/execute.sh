#!/bin/sh
# A bash script to automate the task of running simulations for the segregation of an architecture.
numberOfMonomers=$1
architecture=$2
numberOfRuns=50
runIndex=$3 # The index of the run to be launched;
special_simulation=$4 # Argument for running the simulations in a special simulation (initialization procedure) folder

# Read before use: Pre-requisites before running this script:
# - The initial mixed states for a particular architecture must have been generated already
# - A folder named as the architecture of interest should exist under b200/ or b500/ and must contain an execute.sh script that can start different runs with different random seeds. For an example, see b200/Arc2-2/execute.sh

# Getting path to base folder:
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd $SCRIPT_DIR/../../ && pwd)/"

SEGREGATION="${BASE_DIR}Segregation/"
SCRIPTS="${BASE_DIR}Global_Scripts/"

# Checking the inputs:
if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] || [ -z ${runIndex} ] || [ -z ${special_simulation} ] # -z flag checks if the variable has zero length
then
	echo "Not enough arguments passed! Please pass the number of monomers, the architecture, the run index, and the name of the initialization procedure as arguments while invoking the script"
	echo "For example: $ execute.sh 200 Arc-0 1 fene_recenter"
	exit 1
fi

#echo ${architecture}
# NOTE: Have the execute.sh ready with the random seeds for this architecture!
# Making the necessary directories:
cd ${SEGREGATION}
if ! [ -z ${special_simulation} ]; then
	folderSuffix=b${numberOfMonomers}/${special_simulation}/
	echo "Copying files to the special simulation folder '${folderSuffix}'."
else
	folderSuffix=b${numberOfMonomers}/
	echo "No special simulation folder passed. Copying files in the default Segregation/b${numberOfMonomers} directory."
fi

# navigating to the relevant directory:
error=0
cd ${folderSuffix}
mkdir -p ${architecture}
error=$[${error}+$?] # adding the exit code of the previous command
# Making the run directories:
#for((i=1; i<=numberOfRuns; i++)) # Can be done for multiple runs if desired
mkdir -p ${architecture}/run${runIndex}
cd ${SEGREGATION}Scripts/

# copying the initial mixed states that should have been generated already:
echo "Copying mixed states files to Segregation/${folderSuffix}."
# Makign directory where the mixed states will be copied:
mkdir -p ${architecture}/run${runIndex}/
gcc copyMixedState.c -lm -std=gnu99 -o copy.out -I $BASE_DIR
./copy.out ${numberOfMonomers} ${architecture} 0 ${special_simulation} ${SEGREGATION} ${runIndex} # The 0 indicates that the polymers should not be separated into different files; the copyMixedState.c script will copy the mixed state file as it is to the relevant run folder
error=$?
error=0
if ! [ ${error} -eq 0 ]
then
		echo "Something went wrong while copying the mixed states. Exit code: ${error}"
		exit ${error}
fi
echo "Finished copying mixed states for ${architecture}"
error=0
# running the simulation? Will need to wait for it to get over
# copying LAMMPS script:
#cd b${numberOfMonomers}/
#error=$[${error}+$?]
#cp segregate.lammps ${architecture}/
#error=$[${error}+$?]
#cd ${architecture}
#error=$[${error}+$?]


if ! [ ${error} -eq 0 ]
then
		echo "Something went wrong while copying the segregation.lammps script. Exit code: ${error}"
		exit ${error}
fi

# # copying the execute script in case a special simulation is being used:
# if ! [ -z ${special_simulation} ]; then
# 	cp ${SEGREGATION}b${numberOfMonomers}/${architecture}/execute.sh ../${folderSuffix}${architecture}/
# 	error=$?
# 	if ! [ ${error} -eq 0 ]; then
# 		echo "Something went wrong while copying the execute.sh file to the ${folderSuffix} folder!"
# 		exit ${error}
# 	fi
# 	echo "Copied the execute.sh script to the corresponding ${folderSuffix} folder."
# fi

# Reading the confinement diameter for the current architecture:
. ${SCRIPTS}bash_functions.sh # 'Importing' the bash functions file
diameterFileDir=${SEGREGATION}b${numberOfMonomers}/
if [[ ${special_simulation} == *"inf"* ]]; then # If the segregation needs to be done in an infinite cylinder
	diametersFile=${diameterFileDir}c_Diameters.csv # Path to the diameter database file for the infinite cylinder (constant confinement) simulations
else
	diametersFile=${diameterFileDir}Diameters.csv # Path to the diameter database file for the finite cylinder simulations
fi
 # An empty dictionary that will be filled:
declare -A diameters=()
read_diameters ${diametersFile} # Calling the function from bash_functions
radius=$(bc <<< "scale=2; ${diameters[${architecture}]}/2")
echo "Confinement radius for ${architecture} = ${radius}"
error=$?
if ! [ ${error} -eq 0 ] # checking if the above process ended in an error
then
        echo "A problem was encountered while reading the confinement diameters! Terminating."
        exit ${error}
fi

# Moving to the architecture directory and giving runs:
cd ${SEGREGATION}${folderSuffix}/${architecture}/

# Setting the path to the LAMMPS script:
SEG_SCRIPT=${BASE_DIR}Segregation/Scripts/fene_segregate.lammps
if [[ ${special_simulation} == *"inf"* ]]; then # If the segregation needs to be done in an infinite cylinder
	SEG_SCRIPT=${BASE_DIR}Segregation/Scripts/inf_segregate.lammps
fi

# Defining random seeds for the runs:
seeds=(190686788 260874575 747983061 142559277 635050179 582691149 149585093 820715049 693014654 591232730 396476315 853559767 682319630 221122296 96423063 435689196 334467652 578248473 575945357 558069167 83642928 707654260 688590199 76985948 890405116 263684786 673465517 839630143 515793710 757108445 732602762 476633237 683112204 321051233 555273437 769601417 798673965 282100765 766468201 239825308 755619018 49558515 782161277 308216577 732353746 629267811 597281482 468121962 480789597 357648192)

# LAMMPS_EXEC=`which lmp`
LAMMPS_EXEC=`which lmp`
if [ -z ${LAMMPS_EXEC} ]; then
	echo "LAMMPS executable lmp not found. Please check if the executable lmp is in your PATH, or if the executable is named differently."
	exit 1
fi
i=$runIndex
# Multiple runs can be run simultaneously if desired:
#for((i=1; i<=50; i++))
#do
cd run$i
# using nohup to launch the simulation processes in the background
nohup ${LAMMPS_EXEC} -in ${SEG_SCRIPT} -var seed ${seeds[$[$i-1]]} -var numberOfMonomers ${numberOfMonomers} -var radius ${radius} -log log.lammps > nohup.out &
sleep 1
# done
