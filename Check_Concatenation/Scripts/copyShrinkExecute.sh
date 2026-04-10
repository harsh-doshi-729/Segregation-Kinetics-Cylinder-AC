#!/bin/sh
# This script copies the final mixed state that is produced after a shrink_relax run in the Create_Initial_States/b200/ or Create_Initial_States/b500/ directory.
# It then runs the segregate script(s) to check for concatenations.
# In summary, the LAMMPS script simulates the polymers in a big box for a long time. If the polymers drift apart, they are unconcatenated.
# Similarly, for architectures with multiple small loops, the script simulates the loops only and checks if they drift apart.

# Accepting command line arguments:
numberOfMonomers=$1
architecture=$2
runIndex=$4
flag=$3 # A flag to indicate whether the polymers should be checked for internal concatenations within each polymer
procedure=$5 # Argument for specifying the initialization procedure that was used for creating the mixed state.
# This is the name of the subdirectory where the initialization simulation was run

# Integers to seed the random number generators during the simulations (50 seeds total):
seeds=(190686788 260874575 747983061 142559277 635050179 582691149 149585093 820715049 693014654 591232730 396476315 853559767 682319630 221122296 96423063 435689196 334467652 578248473 575945357 558069167 83642928 707654260 688590199 76985948 890405116 263684786 673465517 839630143 515793710 757108445 732602762 476633237 683112204 321051233 555273437 769601417 798673965 282100765 766468201 239825308 755619018 49558515 782161277 308216577 732353746 629267811 597281482 468121962 480789597 357648192)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd $SCRIPT_DIR/../../ && pwd)/"

SEGREGATION=${BASE_DIR}Segregation/
CHECK_CONCATENATION=${BASE_DIR}Check_Concatenation/

# Checking if input was given
if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] || [ -z ${flag} ] || [ -z ${runIndex} ] || [ -z ${procedure} ]; then
	echo "Not enough arguments passed! Please pass the number of monomers, the architecture," \
	"a flag to indicate whether the polymers should be checked for internal concatenations, the run number for the initialization run," \
	"and name of the initialization procedure used."
	echo "For example: $ copyShrinkExecute.sh 200 Arc-1-2 1 1 fene_recenter"
	exit 1
fi

# copying file using the copyMixedState.c script
cd ${CHECK_CONCATENATION}
if ! [ -z ${procedure} ]; then
	mkdir -p b${numberOfMonomers}/${procedure}/${architecture}/run${runIndex}
	folderPrefix=b${numberOfMonomers}/${procedure}/ # The folder where the architecture directories reside
	echo "Copying files to the special simulation folder '${folderPrefix}'."

	localPath=${CHECK_CONCATENATION}
fi

# Compiling and running the copyMixedState.c script to copy the mixed state file to the relevant subdirectory in the Check_Concatenation directory:
cd ${SEGREGATION}Scripts/
gcc copyMixedState.c -lm -std=gnu99 -o ${CHECK_CONCATENATION}Scripts/copy.out -I $BASE_DIR
cd ${CHECK_CONCATENATION}Scripts/
if [ ${flag} -eq 1 ]; then # If the internal concantenations of each polymer are to be checked:
	# copying each polymer from the mixed state to a separate file:
	./copy.out ${numberOfMonomers} ${architecture} ${flag} ${procedure} ${CHECK_CONCATENATION} ${runIndex}
fi
# Copying the mixed state without separating the polymers:
./copy.out ${numberOfMonomers} ${architecture} 0 ${procedure} ${CHECK_CONCATENATION} ${runIndex}

# Checking if the copy failed:
exitCode=$?
if ! [ ${exitCode} -eq 0 ]; then
	echo "Copying failed. Terminating."
	exit 1
fi

# Running the segregation simulation scripts:
cd ${CHECK_CONCATENATION}b${numberOfMonomers}/${procedure}/${architecture}/run${runIndex}/
# LAMMPS_EXEC=`which lmp` # Assuming lmp is the LAMMPS executable; change if it is different
LAMMPS_EXEC=`which lmp`
if [ -z ${LAMMPS_EXEC} ]; then
	echo "LAMMPS executable lmp not found. Please check if the executable lmp is in your PATH, or if the executable is named differently."
	exit 1
fi
# Using nohup to launch the simulation processes in the background
if [ ${flag} -eq 1 ]; then # If the internal concantenations of each polymer are to be checked:
	echo "Checking concatenations between different side loops of the same polymer of ${architecture} Run ${runIndex}."
	for((n=1; n<=2; n++)) # Looping over the two polymers in the system; checking internal concantenations individually
	do
		nohup ${LAMMPS_EXEC} -in ${CHECK_CONCATENATION}Scripts/separate_segregate.lammps -var numberOfMonomers ${numberOfMonomers} \
		 -var polymerIndex ${n} -var seed ${seeds[$[$i-1]]} -log separate_log.lammps > "nohup${n}.out" &
		sleep 1
	done
fi
echo "Checking concatenation between two polymers of ${architecture} Run ${runIndex}."    
nohup ${LAMMPS_EXEC} -in ${CHECK_CONCATENATION}Scripts/segregate.lammps -var numberOfMonomers ${numberOfMonomers} \
	-var seed ${seeds[$[$i-1]]} -log both_log.lammps > nohup.out &
sleep 1
