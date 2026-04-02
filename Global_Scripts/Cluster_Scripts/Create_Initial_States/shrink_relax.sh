#!/bin/bash
# This script runs the the shrink-relax + fixing-bonding algorithm given the appropriate architecture and number of monomers. This process is followed by calculating the monomer distributions in different sections and the bond energy distributions.
# Run this script in the background using nohup as done in the execute_shrink_relax.sh script for automation.

numberOfMonomers=$1
architecture=$2
run=$3
# NOTE: Please ensure that the Harsh/Scripts/System_File_Paths/system_file_paths.h config file has the appropriate special simulation set
# orientation angles for each polymer (in degrees):
angle1=$4
angle2=$5
seed=$6
special_simulation=$7 # Optional argument for running the algorithms in a special folder
bond_coeff=1000

# filePaths:
SCRIPTS=/scratch/Harsh/Scripts/
NEW_SEGREGATION=/scratch/Harsh/New_Segregation/
CREATE_INITIAL_STATES=${NEW_SEGREGATION}Create_Initial_States/

if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] || [ -z ${run} ] || [ -z ${angle1} ] || [ -z ${angle2} ] || [ -z ${seed} ] # -z flag checks whether the variable has zero length
then
        echo "Not enough arguments passed! Please pass number of monomers, the architecture, the run number, the two oritentation angles, and the seed  while executing this script!"
        echo "An optional argument for the special simulation identifier can be passed as well to store the files in a special directory"
        exit 1
fi

# A dictionary of seeds where the architecture are the keys:
declare -A seeds=(["Arc0"]=90907660 ["Arc2"]=699158960 ["Arc2-2"]=222493186 ["Arc3"]=108634382 ["Arc4"]=427365698 ["Arc5"]=503472353 ["Arc6"]=609131881 ["Arc7"]=86300493 ["Arc8"]=700039541 ["Arc9"]=785459453 ["Arc10"]=683299087 ["Arc11"]=454081557 ["Arc1_10"]=450923791 ["Arc1_1"]=432295288 ["Arc_Shree_2_5"]=306495922 ["ArcI-8"]=303808738 ["Arc1"]=303791268)
# Note: Mistakenly, Arc10 used the same seed as Arc1
# extra seeds: 303791268 236790580 86793632 675398236

# A dictionary to store axis lengths of the confining cylinder for various architectures:
declare -A diameters=()
# CSV File containing axis lengths:
diametersFile="${NEW_SEGREGATION}b${numberOfMonomers}/Diameters.csv"


cd ${CREATE_INITIAL_STATES}
# Making the directory:
if ! [ -z ${special_simulation} ]; then
	mkdir -p b${numberOfMonomers}/Previous_Attempts/${special_simulation}/${architecture}/run${run}
	folderPrefix=b${numberOfMonomers}/Previous_Attempts/${special_simulation}/ # The folder where the architecure directories reside
	shrink_relax_script=${CREATE_INITIAL_STATES}${folderPrefix}shrink_relax.lammps # The path to the shrink_relax LAMMPS script 
	echo "shrink_relax script: ${shrink_relax_script}"
	echo "Copying files to the special simulation folder 'b${numberOfMonomers}/Previous_Attempts/${special_simulation}'."
else
	mkdir -p b${numberOfMonomers}/${architecture}/run${run}
	folderPrefix=b${numberOfMonomers}/
	shrink_relax_script=${CREATE_INITIAL_STATES}Scripts/fene_shrink_relax.lammps
	# if [ "${architecture}" = "Arc_Lin" ]; then
	# 	shrink_relax_script=${CREATE_INITIAL_STATES}Scripts/linear_mixing.lammps # This does not perfrom shrink-relax, but initializes the linear polymer mixed state
	# fi
	echo "No special simulation folder passed. Copying files in the default Create_Intial_States/b${numberOfMonomers} directory."
fi

# Creating the initial file:
# if [ "${architecture}" != "Arc_Lin" ]; then # skipping for linear architecture
cd Scripts/
gcc createInitialState.c -lm -std=gnu99
./a.out ${numberOfMonomers} ${architecture} ${run}
error=$?
if ! [ ${error} -eq 0 ] # checking if the above process ended in an error
then
		echo "The initial configuration could not be copied successfully! Terminating."
		exit ${error}
fi
# else
# 	echo "Skipping creating initial config file for Arc_Lin. Assuming the config file is already uploaded."
# fi

# Compiling the C scripts for analysis:
runFolder=${CREATE_INITIAL_STATES}${folderPrefix}${architecture}/run${run}/ # storing the compiled executables separately for each run
# Monomer distribution:
cd ${NEW_SEGREGATION}Scripts/
gcc sectional_monomer_density.c -lm -std=gnu99 -o ${runFolder}sect.out
echo "The regions are read from the regions_config.h file in Harsh/Scripts/Config_Files/. Please ensure this information is correctly set!"
gcc general_monomer_density.c -lm -std=gnu99 -o ${runFolder}ge.out
gcc radial_monomer_density.c -lm -std=gnu99 -o ${runFolder}radial.out
# CoM Time series:
gcc calculateCoMTimeSeries.c -lm -std=gnu99 -o ${runFolder}com.out
# Bond distribution: not useful for FENE bonds
cd ${CREATE_INITIAL_STATES}Scripts/
# gcc calculateBondEnergyDistribution.c -lm -std=gnu99 -o ${runFolder}energy.out
# Pair correlation:
gcc pair_correlation.c -lm -std=gnu99 -o ${runFolder}pair.out

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
module load codes/gcc/lammps/2018.03.16
cd run${run}/
echo "Running shrink-relax LAMMPS script"
# NOTE: This script does not run the simulation in the background. It is expected that this entire script will run in the background
echo "Seed = ${seed}"

# Running the shrink-relax simulation:
LAMMPS_EXEC=`which lmp_gcc`
ln -sf ${LAMMPS_EXEC} ${CREATE_INITIAL_STATES}Scripts/Harsh_shrink-relax # Making a symbolic link; pseudoname
#lmp_gcc -in ../../shrink_relax.lammps -var seed ${seeds[$[$i-1]]} -var run $(printf %02d $i) -log temp_log.lammps
# lmp_gcc -in ${CREATE_INITIAL_STATES}Scripts/shrink_relax.lammps -var seed ${seed} -var run $(printf %02d ${run}) -var angle1 ${angle1} -var angle2 ${angle2} -log log.lammps
#lmp_gcc -in ../../cylinder_mixing_200.lammps -var seed ${seeds[$[$i-1]]} -var run $(printf %02d $i)
${CREATE_INITIAL_STATES}Scripts/Harsh_shrink-relax -in ${shrink_relax_script} -var seed ${seed} -var numberOfMonomers ${numberOfMonomers} -var run $(printf %02d ${run}) -var angle1 ${angle1} -var angle2 ${angle2} -var finalRadius ${final_radius} -log log.lammps
error=$?
#echo "The shrink-relax LAMMPS script ended with an error code ${error}"
if ! [ ${error} -eq 0 ] # checking if the above process ended in an error
then
        echo "The shrink-relax LAMMPS script did not run successfully! Terminating."
        exit ${error}
fi

# Running the short simulation to equilibriate the state
# lmp_gcc -in ../../../Scripts/equilibriate.lammps -var seed ${seeds[${architecture}]} -var numberOfMonomers ${numberOfMonomers} -log dist_log.lammps
error=$?
if ! [ ${error} -eq 0 ] # checking if the above process ended in an error
then
        echo "The equilibriate LAMMPS script did not run successfully! Terminating."
        exit ${error}
fi

# Continuing to the distribution analysis:
echo ""
echo "Monomer Distribution:"
mkdir -p monomer_distribution
mkdir -p bond_energy
mkdir -p pair_correlation
# Monomer distribution:
cd ${NEW_SEGREGATION}Scripts/ # going to New_Segregation/Scripts
# Ensure that the Section.h file has the desired sections

binWidth=0.2
${runFolder}sect.out ${numberOfMonomers} ${architecture} ${run} ${binWidth}
rm ${runFolder}sect.out # removing the executable

# calculating regionwise monomer distribution: (regions of different polymer)
echo ""
echo "Regional Monomer Distribution:"
binWidth=0.1 # updating binWidth for radial distribution
axisLength=$( bc <<< "scale=2; 10 * ${final_radius}")
${runFolder}ge.out ${numberOfMonomers} ${architecture} ${run} Create_Initial_States 0 ${axisLength}
rm ${runFolder}ge.out

# calculating radial distribution:
echo ""
echo "Radial Monomer Distribution:"
${runFolder}radial.out ${numberOfMonomers} ${architecture} ${run} ${binWidth} Create_Initial_States
rm ${runFolder}radial.out

# Region wise CoM Time series:
echo ""
echo "Region wise CoM Time Series:"
${runFolder}com.out ${numberOfMonomers} ${architecture} ${run} Create_Initial_States
rm ${runFolder}com.out

# calculating the pair_correlation functions:
echo ""
echo "Pair Correlation Function:"
range=$( ceil ${final_radius} ) # ceil() function imported from bash functions; range is maximum distance up to which the pair_correlation is to be calculated
${runFolder}pair.out ${numberOfMonomers} ${architecture} ${run} ${range} ${binWidth} 1 1
# The second last argument 1 is a flag to indicate that the pair correlations should be calculated separately
# The last argument 1 is a flag to indicate that the bonded neighbours should be ignored while calculating the pair correlation components
rm ${runFolder}pair.out

# Bond Energy Distribution: Not useful for FENE bonds yet
#echo ""
#echo "Bond Energy Distribution:"
#cd ../Create_Initial_States/Scripts/ # going to Create_Initial_States/Scripts/
#${runFolder}energy.out ${numberOfMonomers} ${architecture} ${run} ${bond_coeff}
#rm ${runFolder}energy.out

echo "Done"
