#!/bin/sh
# A bash script to automate the task of running simulations for the segregation of an architecture.
architecture=$2
numberOfRuns=50
numberOfMonomers=$1
SEGREGATION="/Segregation/"
SCRIPTS="/Global_Scripts/"
special_simulation=$3 # Argument for running the simulations in a special simulation (initialization procedure) folder

# Read before use: Pre-requisites before running this script:
# - The initial mixed states for a particular architecture must have been generated already
# - A folder named as the architecture of interest should exist under b200/ or b500/ and must contain an execute.sh script that can start different runs with different random seeds. For an example, see b200/Arc2-2/execute.sh
# - If the single_state special simulations are to be run, then the setup_single_simulation.sh script should be run first

# Checking the inputs:
if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] # -z flag checks if the variable has zero length
then
	echo "Not enough arguments passed! Please pass the number of monomers and the architecture as arguments while invoking the script"
	echo "For example: $ execute.sh 200 Arc0"
	echo "An optional argument for the special simulation identifier can be passed as well to store the files in a special directory"
	exit 1
fi

#echo ${architecture}
# NOTE: Have the execute.sh ready with the random seeds for this architecture!
# Making the necessary directories:
cd ${NEW_SEGREGATION}
if ! [ -z ${special_simulation} ]; then
	# mkdir -p b${numberOfMonomers}/Previous_Attempts/${special_simulation}/ # The fodler should already exist with the segregate.lammps file
	folderSuffix=b${numberOfMonomers}/Previous_Attempts/${special_simulation}/
	echo "Copying files to the special simulation folder '${folderSuffix}'."
else
	folderSuffix=b${numberOfMonomers}/
	echo "No special simulation folder passed. Copying files in the default New_Segregation/b${numberOfMonomers} directory."
fi

if [ "${special_simulation}" != "single_state_low" ] && [ "${special_simulation}" != "single_state_high" ]; then
# When this condition is satisfied, it is the usual case of segregating different mixed states
	# navigating to the relevant directory:
	error=0
	cd ${folderSuffix}
	mkdir -p ${architecture}
	error=$[${error}+$?] # adding the exit code of the previous command
	# creating the folders for the architecture:
	if ! [ "${architecture}" = "Arc0" ] || ! [ -z ${special_simulation} ] # as long as it isnt the default Arc0 folder
	then
		cp ${NEW_SEGREGATION}b${numberOfMonomers}/Arc0/copy.sh ${architecture}/ # copying
	fi
	error=$[${error}+$?] # adding the exit code of the previous command
	cd ${architecture}/
	error=$[${error}+$?] # adding the exit code of the previous command
	bash copy.sh
	error=$[${error}+$?] # adding the exit code of the previous command

	if ! [ ${error} -eq 0 ]
	then
		echo "Something went wrong while creating the runs folders. Exit code: ${error}"
		exit ${error}
	fi

	error=0

	# assuming the execute file has been set up with the different seeds
	cd ${NEW_SEGREGATION}Scripts/

	# copying the initial mixed states that should have been generated already:
	echo "Copying mixed states files to New_Segregation/${folderSuffix}."
	gcc copyMixedState.c -lm -std=gnu99 -o copy.out
	./copy.out ${numberOfMonomers} ${architecture} 0 ${NEW_SEGREGATION}
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

	# copying the execute script in case a special simulation is being used:
	if ! [ -z ${special_simulation} ]; then
		cp ${NEW_SEGREGATION}b${numberOfMonomers}/${architecture}/execute.sh ../${folderSuffix}${architecture}/
		error=$?
		if ! [ ${error} -eq 0 ]; then
			echo "Something went wrong while copying the execute.sh file to the ${folderSuffix} folder!"
			exit ${error}
		fi
		echo "Copied the execute.sh script to the corresponding ${folderSuffix} folder."
	fi
else # This is when only a single mixed state is desired to be segregated in independent runs
	echo "Skipping copying the initial state files and execute file using the normal way. Use setup_single_simulation.sh instead."
fi

# Reading the confinement diameter for the current architecture:
. ${SCRIPTS}bash_functions.sh # 'Importing' the bash functions file
diametersFile=${NEW_SEGREGATION}b${numberOfMonomers}/Diameters.csv # Path to the diameter database file
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
cd ${NEW_SEGREGATION}${folderSuffix}/${architecture}/
bash execute.sh ${radius}

error=$?
if ! [ ${error} -eq 0 ]; then
	echo "Something went wrong while running the execute.sh files!"
fi
