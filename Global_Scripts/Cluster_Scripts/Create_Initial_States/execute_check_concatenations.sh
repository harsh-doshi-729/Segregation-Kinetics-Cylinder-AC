#!/bin/bash
# This script copies the corresponding mixed states generated in this folder to the Check_Concatenation folder and runs the simulation to check for concatenations. It does this by invoking the script in the Check_Concatentation/Scripts/ folder.
# In order to not overwhelm the system with lots of runs, the script uses the wait command to execute new runs only after a set of them has finished. During "waiting", the terminal is suspended.
# It is advisable to run this script with nohup.

numberOfMonomers=$1
architecture=$2
special_simulation=$3 # Optional argument for running the algorithms in a special folder
numberOfRuns=50
if [ -z ${numberOfMonomers} ] || [ -z ${architecture} ] # -z flag checks whether the variable has zero length
then
        echo "Not enough arguments passed! Please pass number of monomers and the architecture while executing this script!"
		echo "An optional argument for the special simulation identifier can be passed as well to store the files in a special directory"
        exit 1
fi

cd /scratch/Harsh/New_Segregation/Check_Concatenation/Scripts/

internalLoopsFlag=0 # A flag to indicate whether concatenations between internal loops should be checked

# checking for concatenation between two polymers:
PIDs=()
for((i=1; i<=${numberOfRuns}; i++))
do
	echo "Run $i"
	bash copyShrinkExecute.sh ${numberOfMonomers} ${architecture} 0 $i ${special_simulation} & # running for state 1 and polymers unseparated
	sleep 1
	PID="$!"
	if [ -z ${PID} ]; then
		echo "PID not set"
		exit 2
	fi
	echo "PID: ${PID}"
	PIDs+=($PID) # PID of most recent background process
done
# bash copyShrinkExecute.sh ${numberOfMonomers} ${architecture} 0 50 ${special_simulation} # running for state 50 and polymers unseparated
echo "PIDs: ${PIDs[@]}"
wait ${PIDs[@]} # waiting for all those processes to finish so that the cluster is not overwhelmed with tasks

# checking for concatenation within internal loops:
# architectures with more than one internal loop:
internalLoopArchs=("Arc2" "Arc2-2" "Arc6" "Arc7" "Arc8" "Arc9" "Arc10" "Arc11" "Arc1_10" "Arc_Shree_2_5" "Arc_Loop_5")

PIDs=()
if [[ "${internalLoopArchs[@]}" =~ "${architecture}" ]] # Checking if the architecture exists in the array
then
	internalLoopsFlag=1
	for((i=1; i<=${numberOfRuns}; i++))
	do
		echo "Run $i - Internal concatenation"
		bash copyShrinkExecute.sh ${numberOfMonomers} ${architecture} 1 $i ${special_simulation} & # running for state 1 and polymers separated
	sleep 1
	PID=$!
	PIDs+=("${PID}")
	done
	#bash copyShrinkExecute.sh ${numberOfMonomers} ${architecture} 1 50 ${special_simulation} # running for state 50 and polymers separated
fi

wait ${PIDs[@]}
# predicting concatenations across polymers:
bash predict_concatenation.sh ${numberOfMonomers} ${architecture} ${internalLoopsFlag} ${special_simulation}

