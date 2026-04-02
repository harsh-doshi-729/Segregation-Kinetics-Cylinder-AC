#!/bin/sh
# This script copies the final mixed state that is produced after a shrink_relax run in the Create_Initial_States/b200/ or Create_Initial_States/b500/ directory. It then runs the segregate script to check for concatenations.

numberOfMonomers=$1
architecture=$2
runIndex=$4
flag=$3
special_simulation=$5 # Optional argument for running the algorithms in a special folder
seeds=(190686788 260874575 747983061 142559277 635050179 582691149 149585093 820715049 693014654 591232730 396476315 853559767 682319630 221122296 96423063 435689196 334467652 578248473 575945357 558069167 83642928 707654260 688590199 76985948 890405116 263684786 673465517 839630143 515793710 757108445 732602762 476633237 683112204 321051233 555273437 769601417 798673965 282100765 766468201 239825308 755619018 49558515 782161277 308216577 732353746 629267811 597281482 468121962 480789597 357648192)

NEW_SEGREGATION=/scratch/Harsh/New_Segregation/
CHECK_CONCATENATION=${NEW_SEGREGATION}Check_Concatenation/

# Checking if input was given
if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] || [ -z ${flag} ] || [ -z ${runIndex} ]; then
	echo "Not enough arguments passed! Please pass the number of monomers, the architecture, a flag to indicate whether the polymers should be separated to different files, and the run number for the shrink_relax run."
	echo "An optional argument for the special simulation identifier can be passed as well to store the files in a special directory"
	exit 1
fi

# copying file using the copyMixedState.c script
cd ${CHECK_CONCATENATION}
if ! [ -z ${special_simulation} ]; then
	mkdir -p b${numberOfMonomers}/Previous_Attempts/${special_simulation}/${architecture}/run${runIndex}
        folderPrefix=b${numberOfMonomers}/Previous_Attempts/${special_simulation}/ # The folder where the architecture directories reside
        echo "Copying files to the special simulation folder '${folderPrefix}'."

	localPath=${CHECK_CONCATENATION}
else
	echo "No special simulation folder passed. Copying files in the default Check_Concatenation/ directory."
	mkdir -p b${numberOfMonomers}/${architecture}/run${runIndex}/
	localPath=${CHECK_CONCATENATION}
fi

# navigating to New_Segregation/Scripts
cd ${NEW_SEGREGATION}Scripts/
gcc copyMixedState.c -lm -std=gnu99 -o ../Check_Concatenation/Scripts/copy.out
cd ../Check_Concatenation/Scripts/
./copy.out ${numberOfMonomers} ${architecture} ${flag} ${CHECK_CONCATENATION} ${runIndex}
# cp ../Create_Initial_States/b${numberOfMonomers}/${architecture}/langevin_mixed_state_$(printf %02d ${runIndex}).txt b${numberOfMonomers}/${architecture}/run${runIndex}/initial_configuration.txt

# Checking if the copy failed:
exitCode=$?
if ! [ ${exitCode} -eq 0 ]; then
	echo "Copying failed. Terminating."
	exit 1
fi

# executing the segregate scripts:
bash execute.sh ${numberOfMonomers} ${architecture} ${flag} ${runIndex} ${special_simulation}
# cd b${numberOfMonomers}/${architecture}/run${runIndex}/
# module load codes/gcc/lammps/2018.03.16
# nohup lmp_gcc -in ../../../segregate.lammps -var numberOfMonomers ${numberOfMonomers} -var seed ${seeds[$[$i-1]]} > nohup.out &
# sleep 1

