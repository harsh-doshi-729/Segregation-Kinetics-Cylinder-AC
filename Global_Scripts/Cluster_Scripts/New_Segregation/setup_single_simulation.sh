#!/bin/bash
# This script is used to set up the files for single_state low and high simulations. 
# A single mixed state from a previously generated set of states is copied and used to give 50 independent segregation runs from this same state.

NEW_SEGREGATION=/scratch/Harsh/New_Segregation/
numberOfRuns=50
# Variables:
numberOfMonomers=$1
architecture=$2
simulation_type=$3 # The name of the simulation type / algorithm with which the desired mixed state was generated
low_runIndex=$4 # The runIndex for a run that showed low segregation time in the original set of simulations
high_runIndex=$5 # The runIndex for a run that showed high segregation time in the original set of simulations

# Checking if all arguments are passed:
if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] || [ -z ${simulation_type} ] || [ -z ${low_runIndex} ] || [ -z ${high_runIndex} ]
then
    echo "Not enough arguments passed! Please pass the following arguments on the command line in order:"
    echo "Number of monomers"
    echo "Architecture name"
    echo "Simulation type: name of simulation used to generate desired mixed state"
    echo "Run Index for a run that showed low segregation time"
    echo "Run Index for a run that showed a high segregation time"
    exit 1
fi


runs=(${low_runIndex} ${high_runIndex})
special_simulation_folders=("single_state_low" "single_state_high")

# Copying low and high states:
for((n=0; n<2; n++))
do
    cd ${NEW_SEGREGATION}
    destination=${special_simulation_folders[$n]}/b${numberOfMonomers}/${architecture}/
    mkdir -p ${destination}
    cp b${numberOfMonomers}/${architecture}/copy.sh ${destination}

    error=$?
    if ! [ ${error} -eq 0 ]; then
        echo "Something went wrong while copying the copy.sh script. Terminating."
        exit ${error}
    fi

    cd ${destination}
    bash copy.sh # Creating the run folders
    cd ${NEW_SEGREGATION}

    source=b${numberOfMonomers}/Previous_Attempts/${simulation_type}/${architecture}/run${runs[$n]}/initial_configuration.txt
    for((i=1; i<=50; i++))
    do
        cp ${source} ${destination}/run$i/
    done
    error=$?
    if ! [ ${error} -eq 0 ]; then
        echo "Something went wrong while copying the initial state files. Terminating."
        exit ${error}
    fi

    cd ${destination}
    # Copying the execute file:
    cp ../execute.sh ./
    error=$?
    if ! [ ${error} -eq 0 ]; then
        echo "Something went wrong while copying the execute.sh file."
    fi
    # Writing the source to a text file:
    echo "Source initial mixed state: ${NEW_SEGREGATION}${source}" > source.txt
    echo "Copied the initial state file from run${runs[$n]} to ${special_simulation_folders[$n]}."
done

