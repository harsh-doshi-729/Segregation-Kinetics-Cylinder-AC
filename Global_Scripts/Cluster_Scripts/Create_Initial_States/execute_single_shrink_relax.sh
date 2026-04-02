#!/bin/sh
# This script accepts the inputs of the number of monomers and architecture from the command line and performs the shrink-relax+fixing-bonding (fbsr) algorithm as the first step to creating mixed states. This single preliminary mixed state is checked for local density imbalances.
numberOfMonomers=$1
architecture=$2
# NOTE: Please ensure that the Harsh/Scripts/System_File_Paths/system_file_paths.h config file has the appropriate special simulation set

# Orientation angles for the two polymers:
angle1=$3
angle2=$4

special_simulation=$5 # Optional argument for running the algorithms in a special folder
numberOfRuns=1
. /scratch/Harsh/Scripts/bash_functions.sh # Importing script of bash functions
# Moving to the specific directory:
cd /scratch/Harsh/New_Segregation/Create_Initial_States/

if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] # -z flag checks whether the variable has zero length
then
        echo "Not enough arguments passed! Please pass number of monomers, the architecture, and the two orientation angles while executing this script!"
        echo "An optional argument for the special simulation identifier can be passed as well to store the files in a special directory"
        exit 1
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
echo "Seeds: ${seeds[@]}"

# Making the special simulation directory if passed:
if ! [ -z ${special_simulation} ]; then
	mkdir -p b${numberOfMonomers}/Previous_Attempts/${special_simulation}/${architecture}/
        folderPrefix=b${numberOfMonomers}/Previous_Attempts/${special_simulation}/
        echo "Copying files to the special simulation folder 'b${numberOfMonomers}/Previous_Attempts/${special_simulation}'."
else
        mkdir -p b${numberOfMonomers}/${architecture}/
        folderPrefix="b${numberOfMonomers}/"
	echo "No special simulation folder passed. Copying files in the default Create_Initial_States/b${numberOfMonomers} directory."
fi

# Running the shrink-relax script:
cd Scripts/
for((run=1; run<=${numberOfRuns}; run++))
do
	# run=1
	echo "Run${run}"
	mkdir -p ../${folderPrefix}${architecture}/run${run}/
	nohup bash single_shrink_relax.sh ${numberOfMonomers} ${architecture} ${run} ${angle1} ${angle2} ${seeds[$[$run-1]]} ${special_simulation} > ../${folderPrefix}${architecture}/run${run}/nohup.out & # running in the background
#	nohup bash shrink_relax.sh ${numberOfMonomers} ${architecture} ${run} ${angle1} ${angle2} ${seeds[${architecture}]} ${special_simulation} > ../${folderPrefix}b${numberOfMonomers}/${architecture}/run${run}/nohup.out & # running in the background
	sleep 1
done
# i=4

# Old generation of mixed states:
# module load codes/gcc/lammps/2018.03.16
#nohup mpirun -np 4 lmp_gcc -in ../cylinder_mixing_${numberOfMonomers}.lammps -var seed ${seeds[$[$i-1]]} > nohup.out &
#sleep 1
