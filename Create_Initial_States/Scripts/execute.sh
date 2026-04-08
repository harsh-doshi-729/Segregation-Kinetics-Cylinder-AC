#!/bin/bash
# This script runs the shrink-relax procedure starting from two polymers in a big and wide cylinder to create a mixed state of the polymers in a finite cylinder

numberOfMonomers=$1
architecture=$2
run=$3
# NOTE: Please ensure that the Global_Scripts/System_File_Paths/system_file_paths.h config file has the appropriate (desired) special simulation set
seed=$4
special_simulation=$5 # Argument for running a particular initialization algorithm in its special folder

# filePaths:
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd $SCRIPT_DIR/../../ && pwd)/"
SCRIPTS=${BASE_DIR}Global_Scripts/
CREATE_INITIAL_STATES=${BASE_DIR}Create_Initial_States/
SEGREGATION=${BASE_DIR}Segregation/

if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] || [ -z ${run} ] || [ -z ${seed} ] || [ -z ${special_simulation} ] # -z flag checks whether the variable has zero length
then
        echo "Not enough arguments passed! Please pass number of monomers, the architecture, the run number, the seed, and the initialization procedure folder name while executing this script!"
#         echo "An optional argument for the special simulation identifier can be passed as well to store the files in a special directory"
        exit 1
fi

# A dictionary to store axis lengths of the confining cylinder for various architectures:
declare -A diameters=()
# CSV File containing diameters and axis lengths of the cylinder for a chosen polymer architectures:
diametersFile="${SEGREGATION}b${numberOfMonomers}/Diameters.csv"


cd ${CREATE_INITIAL_STATES}
# Making the directory where the simulation will be run; it is named based on the initialization procedure and architecture chosen:
if ! [ -z ${special_simulation} ]; then
	mkdir -p b${numberOfMonomers}/${special_simulation}/${architecture}/run${run}
	folderPrefix=b${numberOfMonomers}/${special_simulation}/ # The folder where the architecure directories reside
	echo "shrink_relax script: ${shrink_relax_script}"
	echo "Copying files to the special simulation folder 'b${numberOfMonomers}/${special_simulation}'."
else # The directory for the will be made directly under b<N>/ instead of in a <special_simulation> (initialization procedure) directory
	mkdir -p b${numberOfMonomers}/${architecture}/run${run}
	folderPrefix=b${numberOfMonomers}/
	if [ "${architecture}" = "Arc_Lin" ]; then
		shrink_relax_script=${CREATE_INITIAL_STATES}Scripts/linear_mixing.lammps # This does not perfrom shrink-relax, but initializes the linear polymer mixed state
	fi
	echo "No special simulation folder passed. Copying files in the default Create_Intial_States/b${numberOfMonomers} directory."
fi

# Creating the initial state from which the simulation will be run using the createInitialState.c script:
cd ${CREATE_INITIAL_STATES}Scripts/
gcc createInitialState.c -lm -std=gnu99 -o create.out -I ${BASE_DIR}
./create.out ${numberOfMonomers} ${architecture} ${special_simulation} ${run}
error=$?
if ! [ ${error} -eq 0 ] # checking if the above process ended in an error
then
        echo "A problem was encountered while creating the initial state! Terminating."
        exit ${error}
fi

# shrink_relax_script=${CREATE_INITIAL_STATES}${folderPrefix}shrink_relax.lammps # The path to the shrink_relax LAMMPS script 
shrink_relax_script=${CREATE_INITIAL_STATES}Scripts/${special_simulation}.lammps # The path to the shrink_relax LAMMPS script 

# Compiling the C scripts for analysis:
runFolder=${CREATE_INITIAL_STATES}${folderPrefix}${architecture}/run${run}/ # storing the compiled executables separately for each run
# Monomer distribution:
cd ${SEGREGATION}Scripts/
echo "The regions are read from the regions_config.h file in Global_Scripts/Config_Files/. Please ensure this information is correctly set!"
gcc general_monomer_density.c -lm -std=gnu99 -o ${runFolder}ge.out -I ${BASE_DIR}
gcc radial_monomer_density.c -lm -std=gnu99 -o ${runFolder}radial.out -I ${BASE_DIR}
# CoM Time series:
gcc calculateCoMTimeSeries.c -lm -std=gnu99 -o ${runFolder}com.out -I ${BASE_DIR}

cd ${CREATE_INITIAL_STATES}Scripts/
gcc sectional_monomer_density.c -lm -std=gnu99 -o ${runFolder}sect.out -I ${BASE_DIR}


# Reading the diameter from Diameters.csv database using bash functions defined in another script
. ${SCRIPTS}bash_functions.sh # 'Importing' the bash functions file
read_diameters ${diametersFile} # Calling the function from bash_functions
final_radius=$(bc <<< "scale=2; ${diameters[${architecture}]}/2") # using base conversion tool to perform floating point calculation
error=$?
if ! [ ${error} -eq 0 ] # checking if the above process ended in an error
then
        echo "A problem was encountered while reading the confinement diameters! Terminating."
        exit ${error}
fi
echo "Confinement radius for ${architecture} = ${final_radius}"

# Running simulation:
cd ../${folderPrefix}${architecture}/
cd run${run}/ # Running the simulation in this folder
echo "Running shrink-relax LAMMPS script"
# NOTE: This script does not run the simulation in the background. It is expected that this entire script will run in the background using nohup
echo "Seed = ${seed}"

# Running the shrink-relax simulation:
# LAMMPS_EXEC=`which lmp` # Replace the name of the LAMMPS executable for your installation if lmp is not found
LAMMPS_EXEC=/usr/bin/lmp
if [ -z ${LAMMPS_EXEC} ]; then
	echo "LAMMPS executable not found in PATH! Please ensure LAMMPS is installed and the executable is in the PATH, or update the LAMMPS_EXEC variable in this script to point to the executable."
	exit 1
fi
# Shrink-Relax:
${LAMMPS_EXEC} -in ${shrink_relax_script} -var seed ${seed} -var numberOfMonomers ${numberOfMonomers} -var run $(printf %02d ${run}) -var angle1 ${angle1} -var angle2 ${angle2} -var finalRadius ${final_radius} -log shrink_log.lammps
error=$?
if ! [ ${error} -eq 0 ] # checking if the above process ended in an error
then
        echo "The shrink-relax LAMMPS script did not run successfully! Terminating."
        exit ${error}
fi

# Continuing to the distribution analysis:
echo ""
echo "Monomer Distribution:"
mkdir -p monomer_distribution
# Monomer distribution:
cd ${SEGREGATION}Scripts/ # going to Segregation/Scripts
# Ensure that the Global_Scripts/Utility/Section.h file has the desired sections

if [[ ${special_simulation} == *"inf"* ]]; then # Checking for the presence of "inf" indicating simulation in an infinite cylinder
	axisLength="" # Empty
else
	axisLength=$( bc <<< "scale=2; 10 * ${final_radius}")
fi

binWidth=0.2
${runFolder}sect.out ${numberOfMonomers} ${architecture} ${run} ${binWidth} ${special_simulation} ${axisLength}
rm ${runFolder}sect.out # removing the executable

# calculating regionwise monomer distribution: (regions of different polymer)
echo ""
echo "Regional Monomer Distribution:"
binWidth=0.1 # updating binWidth for radial distribution
${runFolder}ge.out ${numberOfMonomers} ${architecture} ${run} Create_Initial_States ${special_simulation} 0 ${axisLength}
rm ${runFolder}ge.out

# calculating radial distribution:
echo ""
echo "Radial Monomer Distribution:"
${runFolder}radial.out ${numberOfMonomers} ${architecture} ${run} ${binWidth} ${special_simulation} Create_Initial_States
rm ${runFolder}radial.out

# Region wise CoM Time series:
echo ""
echo "Region wise CoM Time Series:"
${runFolder}com.out ${numberOfMonomers} ${architecture} ${run} ${special_simulation} Create_Initial_States
rm ${runFolder}com.out

echo "Done"
