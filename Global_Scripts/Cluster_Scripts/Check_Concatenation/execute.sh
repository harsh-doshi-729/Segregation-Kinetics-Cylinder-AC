#!/bin/sh
# This script runs the simulation of two polymers in a big box to check if they are concatenated or not
module load codes/gcc/lammps/2018.03.16
seeds=(190686788 260874575 747983061 142559277 635050179 582691149 149585093 820715049 693014654 591232730 396476315 853559767 682319630 221122296 96423063 435689196 334467652 578248473 575945357 558069167 83642928 707654260 688590199 76985948 890405116 263684786 673465517 839630143 515793710 757108445 732602762 476633237 683112204 321051233 555273437 769601417 798673965 282100765 766468201 239825308 755619018 49558515 782161277 308216577 732353746 629267811 597281482 468121962 480789597 357648192)
numberOfMonomers=$1
arc=$2
run=$4
flag=$3
special_simulation=$5 # Optional argument for running the algorithms in a special folder

if [ -z ${arc} ] || [ -z ${numberOfMonomers} ] || [ -z ${run} ] || [ -z ${flag} ]; then
        echo "Not enough arguments passed! Please pass the number of monomers, the architecture, the flag to indicate whether the polymers should be separated, and the run number."
		echo "An optional argument for the special simulation identifier can be passed as well to store the files in a special directory"
        exit 1
fi

CHECK_CONCATENATION=/scratch/Harsh/New_Segregation/Check_Concatenation/
cd ${CHECK_CONCATENATION}
# Making the necessary directories:
if ! [ -z ${special_simulation} ]; then
	mkdir -p b${numberOfMonomers}/Previous_Attempts/${special_simulation}/${arc}
        folderPrefix=b${numberOfMonomers}/Previous_Attempts/${special_simulation}/ # The folder where the architecure directories reside
        echo "Copying files to the special simulation folder 'b${numberOfMonomers}/Previous_Attempts/${special_simulation}'."
else
	mkdir -p b${numberOfMonomers}/${arc}/
	folderPrefix=b${numberOfMonomers}/
	echo "No special simulation folder passed. Copying files in the default Check_Concatenation/ directory."
fi

# making symbolic link to LAMMPS executable:
LAMMPS_EXEC=`which lmp_gcc`
ln -sf ${LAMMPS_EXEC} Scripts/Harsh_concat

cd ${folderPrefix}${arc}/
startRun=1
numberOfRuns=5
customRuns=(9 19 26)
lenRuns=3

PIDs=()
# for((i=${startRun}; i<=${numberOfRuns}; i++))
#for((n=0; n<${lenRuns}; n++))
# do
#    i=${customRuns[$n]}
    i=${run}
    mkdir -p run$i
    cd run$i/
    echo ${seeds[$[$i-1]]}
    if [ ${flag} -eq 1 ]; then
	echo "Checking concatenations between different side loops of the same polymer of ${arc} Run ${run}."
	for((n=1; n<=2; n++))
	do
		nohup ${CHECK_CONCATENATION}Scripts/Harsh_concat -in ${CHECK_CONCATENATION}Scripts/separate_segregate.lammps -var numberOfMonomers ${numberOfMonomers} -var polymerIndex ${n} -var seed ${seeds[$[$i-1]]} -log separate_log.lammps > "nohup${n}.out" # not running this in the background
	PIDs+=("$!")
	sleep 1
	done
	cd ..
    else
	echo "Checking concatenation between two polymers of ${arc} Run ${run}."    
	nohup ${CHECK_CONCATENATION}Scripts/Harsh_concat -in ${CHECK_CONCATENATION}Scripts/segregate.lammps -var numberOfMonomers ${numberOfMonomers} -var seed ${seeds[$[$i-1]]} -log both_log.lammps > nohup.out # not running this in the background; assuming this script is run in the background
	PIDs+=("$!")
	sleep 1
	cd ..
    fi
# done
