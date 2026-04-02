#!/bin/sh
# This is a script to automate running the relaxations runs of single polymers to find their radius of gyration (R_g)
# The script must be run using nohup

# input parameters:
numberOfMonomers=$1
architecture=$2

CREATE_INITIAL_STATES=/scratch/Harsh/New_Segregation/Create_Initial_States/
R_G=/scratch/Harsh/New_Segregation/R_g/

if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] # -z flag checks whether the variable has zero length
then
	echo "Not enough arguments passed! Please pass number of monomers and the architecture while executing this script!"
	exit 1
fi
# else:

cd $R_G

# prepare the Arc initial file:
if [ "${architecture}" != "Arc_Lin" ]; then # skipping for linear architecture
	gcc ${CREATE_INITIAL_STATES}Scripts/createInitialState.c -lm -std=gnu99
	./a.out ${numberOfMonomers} ${architecture} ${R_G}
	error=$?
	if ! [ ${error} -eq 0 ] # checking if the above process ended in an error
	then
        	echo "The initial configuration could not be created successfully! Terminating."
	        exit ${error}
	fi
else
	echo "Skipping creating initial config file for Arc_Lin. Assuming the config file is already uploaded."
fi

# Compiling calculation script
gcc ${R_G}calculate_Rg.c -lm -o Rg.out -std=gnu99

# load module
module load codes/gcc/lammps/2018.03.16
LAMMPS_EXEC=`which lmp_gcc`
ln -sf ${LAMMPS_EXEC} Harsh_Rg # creating a symlink of the lmp executable

# launching mc simulation
cd b${numberOfMonomers}/${architecture}/
if [ "${architecture}" == "Arc0" ] || [ "${architecture}" == "Arc_Lin" ]; then # no cross links
	mc_np=4
else
	mc_np=1
fi
seed=197365393 # same seed for all architectures
mpirun -np ${mc_np} ../../Harsh_Rg -in ../../mc_relax.lammps -var numberOfMonomers ${numberOfMonomers} -var seed ${seed}
error=$?
if ! [ ${error} -eq 0 ]; then
	echo "Something went wrong while running the MC LAMMPS script. Terminating."
	exit $error
fi
# launching Langevin simulation
mpirun -np 4 ../../Harsh_Rg -in ../../relax.lammps -var numberOfMonomers ${numberOfMonomers} -var seed ${seed}
error=$?
if ! [ ${error} -eq 0 ]; then
        echo "Something went wrong while running the Langevin LAMMPS script. Terminating."
        exit $error
fi


# calculating R_g after simulation has finished:
# echo "Currently in `pwd`"
cd ../../
./Rg.out ${numberOfMonomers} ${architecture}
