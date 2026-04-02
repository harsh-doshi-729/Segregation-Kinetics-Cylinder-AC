#!/bin/bash
# This is a script to perform the analysis after the shrink-relax variant simulation has been run.
# This script calculates the regional monomer density, the sectional monomer density, the CoM time series,
# the radial monomer density, and the pair correlation.

numberOfMonomers=$1
architecture=$2
run=$3
# NOTE: Please ensure that the Harsh/Scripts/System_File_Paths/system_file_paths.h config file has the appropriate special simulation set
special_simulation=$4 # Optional argument for running the algorithms in a special folder

# filePaths:
SCRIPTS=/scratch/Harsh/Scripts/
NEW_SEGREGATION=/scratch/Harsh/New_Segregation/
CREATE_INITIAL_STATES=${NEW_SEGREGATION}Create_Initial_States/

if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] || [ -z ${run} ]  # -z flag checks whether the variable has zero length
then
    echo "Not enough arguments passed! Please pass number of monomers, the architecture, and the run number while executing this script!"
    echo "An optional argument for the special simulation identifier can be passed as well to store the files in a special directory"
    exit 1
fi

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
	if [ "${architecture}" = "Arc_Lin" ]; then
		shrink_relax_script=${CREATE_INITIAL_STATES}Scripts/linear_mixing.lammps # This does not perfrom shrink-relax, but initializes the linear polymer mixed state
	fi
	echo "No special simulation folder passed. Copying files in the default Create_Intial_States/b${numberOfMonomers} directory."
fi


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

cd ../${folderPrefix}${architecture}/
cd run${run}/
# Continuing to the distribution analysis:
echo ""
echo "Monomer Distribution:"
mkdir -p monomer_distribution
# mkdir -p bond_energy
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
binWidth=0.1 # updating binWidth for radial distribution and pair correlation
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
