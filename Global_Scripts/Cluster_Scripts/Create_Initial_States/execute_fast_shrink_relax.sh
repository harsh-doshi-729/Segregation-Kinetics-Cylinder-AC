#!/bin/bash
# This script accepts the inputs of the number of monomers and architecture from the command line and performs the shrink-relax algorithm as the first step to creating mixed states. This single preliminary mixed state is checked for local density imbalances.
# This script attempts to speed up the initialization process by making the MC simulation common for all runs
numberOfMonomers=$1
architecture=$2
# NOTE: Please ensure that the Harsh/Scripts/System_File_Paths/system_file_paths.h config file has the appropriate special simulation set

# Orientation angles for the two polymers:
angle1=$3
angle2=$4

special_simulation=$5 # Optional argument for running the algorithms in a special folder
skip_MC=false # Indicates whether the Monte Carlo step should be skipped
only_MC=true # Indicates whether only the Monte Carlo simulation for run1 should be run
numberOfRuns=50
. /scratch/Harsh/Scripts/bash_functions.sh # Importing script of bash functions
# Moving to the specific directory:
NEW_SEGREGATION=/scratch/Harsh/New_Segregation/
CREATE_INITIAL_STATES=${NEW_SEGREGATION}Create_Initial_States/
cd ${CREATE_INITIAL_STATES}

echo "It is advised to run this script with nohup."

if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] # -z flag checks whether the variable has zero length
then
        echo "Not enough arguments passed! Please pass number of monomers, the architecture, and the two orientation angles while executing this script!"
        echo "An optional argument for the special simulation identifier can be passed as well to store the files in a special directory"
        exit 1
fi

if [ "${architecture}" == "Arc_Lin" ]; then
    echo "The architecture Arc_Lin cannot be sped up with this script yet. Please use the usual execute.sh script for creating Arc_Lin mixed states."
    # Arc_Lin does not have a divided mixing LAMMPS script with separate MC and Langevin dynamics modules.
    # So, both modules must be run for all runs independently
fi

if [ -z ${angle1} ] || [ -z ${angle2} ]; then
        echo "No argument(s) were passed for atleast one of the orientation angles! Try again."
        exit 1
fi

# Reading seeds from the seeds database:
read_seeds seeds.txt ${architecture}
error=$?
if ! [ ${error} -eq 0 ]; then
        echo "Something went wrong while reading the seeds from the seeds database."
        exit 1
fi
# echo "Seeds: ${seeds[@]}"

# Making the special simulation directory if passed:
if ! [ -z ${special_simulation} ]; then
    mkdir -p b${numberOfMonomers}/Previous_Attempts/${special_simulation}/${architecture}/
    folderPrefix=b${numberOfMonomers}/Previous_Attempts/${special_simulation}/
    scriptPrefix=${folderPrefix}divided_shrink_relax/
    echo "Copying files to the special simulation folder 'b${numberOfMonomers}/Previous_Attempts/${special_simulation}'."
else
    mkdir -p b${numberOfMonomers}/${architecture}/
    folderPrefix="b${numberOfMonomers}/"
    scriptPrefix=Scripts/divided_shrink_relax/
    echo "No special simulation folder passed. Copying files in the default Create_Initial_States/b${numberOfMonomers} directory."
fi
mc_script=${CREATE_INITIAL_STATES}${scriptPrefix}mc_shrink_relax.lammps # The path to the Monte Carlo shrink_relax LAMMPS script
shrink_relax_script=${CREATE_INITIAL_STATES}${scriptPrefix}shrink_relax.lammps # The path to the shrink_relax LAMMPS script 
long_run_script=${CREATE_INITIAL_STATES}${scriptPrefix}long_run_shrink_relax.lammps # The path to the long run LAMMPS script after shrinking

# Making run folders:
cd ${CREATE_INITIAL_STATES}Scripts/
for((run=1; run<=${numberOfRuns}; run++))
do
	# run=1
	echo "Run${run}"
	mkdir -p ../${folderPrefix}${architecture}/run${run}/
done

# Compiling the C scripts for analysis:
if ! [ "$only_MC" == "true" ]; then
    for((run=1; run<=${numberOfRuns}; run++))
    do
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
    done
fi

# Reading the diameter from Diameters.csv database using bash functions defined in another script
diametersFile="${NEW_SEGREGATION}b${numberOfMonomers}/Diameters.csv"
# empty dictionary:
declare -A diameters=()
read_diameters ${diametersFile} # Calling the function from bash_functions
error=$?
final_radius=$(bc <<< "scale=2; ${diameters[${architecture}]}/2") # using base conversion tool to perform floating point calculation
error=$[$? + $error]
if ! [ ${error} -eq 0 ] # checking if the above process ended in an error
then
        echo "A problem was encountered while reading the confinement diameters! Terminating."
        exit ${error}
fi
echo "Confinement radius for ${architecture} = ${final_radius}"

# Setting up initial_configuration file:
if ! [ "$skip_MC" == "true" ]; then
    gcc createInitialState.c -lm -std=gnu99
    run=1
    ./a.out ${numberOfMonomers} ${architecture} ${run}
    error=$?
    if ! [ ${error} -eq 0 ] # checking if the above process ended in an error
    then
            echo "The initial configuration could not be copied successfully! Terminating."
            exit ${error}
    fi
fi

# Running the MC shrink-relax script for first run:
cd ${CREATE_INITIAL_STATES}${folderPrefix}${architecture}/run1/
module load codes/gcc/lammps/2018.03.16
LAMMPS_EXEC=`which lmp_gcc`
ln -sf ${LAMMPS_EXEC} ${CREATE_INITIAL_STATES}Scripts/Harsh_shrink-relax # Making a symbolic link; pseudoname
# MC simulation:
if [ "${architecture}" == "Arc0" ]; then # no cross links
        mc_np=4
else
        mc_np=1
fi

if ! [ "$skip_MC" == "true" ]; then
    mpirun -np ${mc_np} ${CREATE_INITIAL_STATES}Scripts/Harsh_shrink-relax -in ${mc_script} -var seed ${seeds[0]} -var numberOfMonomers ${numberOfMonomers} -var run 01 -var angle1 ${angle1} -var angle2 ${angle2} -log log.lammps
    error=$?
    if ! [ ${error} -eq 0 ] # checking if the above process ended in an error
    then
            echo "The Monte Carlo shrink-relax LAMMPS script did not run successfully! Terminating."
            exit ${error}
    fi
    if [ "$only_MC" == "true" ]; then
        echo "Done with MC for $architecture Run 1! Terminating."
        exit 0
    fi
fi

# Once the MC simulation is done, copying the MC state for other runs as well
cd ${CREATE_INITIAL_STATES}${folderPrefix}${architecture}/
for((run=2; run<=${numberOfRuns}; run++))
do
    cp run1/mc_configuration.txt run$run/
done
error=$?
if ! [ ${error} -eq 0 ] # checking if the above process ended in an error
then
        echo "Something went wrong while copying the MC configurations. Terminating"
        exit ${error}
fi

# Running the rest of the shrink-relax simulations for all runs:
cd ${CREATE_INITIAL_STATES}Scripts/
for((run=1; run<=${numberOfRuns}; run++))
do
    nohup bash fast_shrink_relax.sh ${numberOfMonomers} ${architecture} ${run} ${seeds[$[$run-1]]} ${special_simulation} > ${CREATE_INITIAL_STATES}${folderPrefix}${architecture}/run${run}/nohup.out & # running in the background
    sleep 1
done
