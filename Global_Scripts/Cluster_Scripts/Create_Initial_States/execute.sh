#!/bin/sh
# DEPRECATED: This script accepts the inputs of the number of monomers and architecture from the command line and runs the necessary steps to create the initial mixed states
numberOfMonomers=$1
architecture=$2
startRun=2
endRun=2
seeds=(190686788 260874575 747983061 142559277 635050179 582691149 149585093 820715049 693014654 591232730 396476315 853559767 682319630 221122296 96423063 435689196 334467652 578248473 575945357 558069167 83642928 707654260 688590199 76985948 890405116 263684786 673465517 839630143 515793710 757108445 732602762 476633237 683112204 321051233 555273437 769601417 798673965 282100765 766468201 239825308 755619018 49558515 782161277 308216577 732353746 629267811 597281482 468121962 480789597 357648192)

if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] # -z flag checks whether the variable has zero length
then
        echo "Not enough arguments passed! Please pass number of monomers and the architecture while executing this script!"
        exit 1
fi
# else:

# setting up the initial file:
mkdir -p b${numberOfMonomers}/${architecture}
gcc createInitialState.c -lm -std=gnu99
./a.out ${numberOfMonomers} ${architecture}
cd ../b${numberOfMonomers}/${architecture}/
# customRuns=(2 3 4 7 8 10 11 13 14 16 18 20 21 25 26 27 28 31 32 33 34 35 36 37 39 40 41 42 43 44 45 46 49 50)
# len=34
module load codes/gcc/lammps/2018.03.16
for((i=${startRun}; i <= ${endRun}; i++))
#for((n=0; n<${len};n++))
do
	#i=${customRuns[$n]}
	# running the simulation:
	mkdir -p run$i
	cd run$i/
	#echo ${seeds[$[$i-1]]}
	#nohup lmp_gcc -in ../../shrink_relax.lammps -var seed ${seeds[$[$i-1]]} -var run $(printf %02d $i) -log temp_log.lammps > temp_nohup.out &
	nohup lmp_gcc -in ../../../Scripts/shrink_relax.lammps -var seed ${seeds[$[$i-1]]} -var run $(printf %02d $i) -log log.lammps > nohup.out &

	#nohup lmp_gcc -in ../../cylinder_mixing_200.lammps -var seed ${seeds[$[$i-1]]} -var run $(printf %02d $i) > nohup.out &
	sleep 1
	cd ..
done
# i=4

# Old generation of mixed states:
# module load codes/gcc/lammps/2018.03.16
#nohup mpirun -np 4 lmp_gcc -in ../cylinder_mixing_${numberOfMonomers}.lammps -var seed ${seeds[$[$i-1]]} > nohup.out &
#sleep 1
