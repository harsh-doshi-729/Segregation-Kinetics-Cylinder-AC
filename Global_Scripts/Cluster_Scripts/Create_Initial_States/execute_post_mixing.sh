#!/bin/bash
# This script runs the extending mixing simulation after the shrink-relax algorithm is done. The independent mixed states are created sequentially at an interval during this simulation.

numberOfMonomers=$1
architecture=$2 
special_simulation=$3 # Optional argument for running the algorithms in a special folder
run=1

if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ]
then
	echo "Not enough arguments passed! Please pass the number of monomers and architecture as arguments to the command line. Terminating."
	echo "An optional argument for the special simulation identifier can be passed as well to store the files in a special directory"
	exit 1
fi

# A dictionary of seeds where the architecture are the keys:
declare -A seeds=(["Arc0"]=90907660 ["Arc2"]=699158960 ["Arc2-2"]=222493186 ["Arc3"]=108634382 ["Arc4"]=427365698 ["Arc5"]=503472353 ["Arc6"]=609131881 ["Arc7"]=86300493 ["Arc8"]=700039541 ["Arc9"]=785459453 ["Arc10"]=683299087 ["Arc11"]=454081557 ["Arc1_10"]=450923791 ["Arc1_1"]=432295288 ["Arc_Shree_2_5"]=306495922 ["ArcI-8"]=303808738 ["Arc1"]=303791268)

CREATE_INITIAL_STATES=/scratch/Harsh/New_Segregation/Create_Initial_States/
cd ${CREATE_INITIAL_STATES}
# Making the necessary directories:
if ! [ -z ${special_simulation} ]; then
	mkdir -p b${numberOfMonomers}/Previous_Attempts/${special_simulation}/${architecture}/run${run}
	cd b${numberOfMonomers}/Previous_Attempts/${special_simulation}/
        folderPrefix=b${numberOfMonomers}/Previous_Attempts/${special_simulation}/ # The folder where the architecure directories reside
        shrink_relax_script=${CREATE_INITIAL_STATES}${folderPrefix}shrink_relax.lammps # The path to the shrink_relax LAMMPS script 
        echo "Copying files to the special simulation folder 'b${numberOfMonomers}/Previous_Attempts/${special_simulation}'."
else
	mkdir -p b${numberOfMonomers}/${architecture}/run${run}
	cd b${numberOfMonomers}/
	echo "No special simulation folder passed. Copying files in the default Create_Initial_State/ directory."
fi

cd ${architecture}/run${run}/

module load codes/gcc/lammps/2018.03.16
# Making a symbolic link:
LAMMPS_EXEC=`which lmp_gcc`
ln -s ${LAMMPS_EXEC} ${CREATE_INITIAL_STATES}Scripts/Harsh_mix
#nohup mpirun -np 4 lmp_gcc -in ${CREATE_INITIAL_STATES}Scripts/cylinder_mixing.lammps -var numberOfMonomers ${numberOfMonomers} -var seed ${seeds[${architecture}]} -log post_log.lammps > nohup_post_mix.out & 
nohup mpirun -np 4 ${CREATE_INITIAL_STATES}Scripts/Harsh_mix -in ${CREATE_INITIAL_STATES}Scripts/fene_cylinder_mixing.lammps -var numberOfMonomers ${numberOfMonomers} -var seed ${seeds[${architecture}]} -var run $(printf %02d ${run}) -log fene_log.lammps > nohup_fene.out &
sleep 1

