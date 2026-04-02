#!/bin/sh
# A script to delete the positions1.dump and positions2.dump file for each run after the monomer desnity is calculated
numberOfMonomers=$1
#architectures=("Arc0" "Arc3" "Arc5" "Arc11" "Arc2" "Arc2-2" "Arc8" "Arc9" "Arc10" "Arc4" "Arc6" "Arc7" "Arc1_10")
#arcArrayLength=13
cd b${numberOfMonomers}/
#for((n=2; n<${arcArrayLength}; n++))
#do
#arc=${architectures[$n]}
arc=$2
cd ${arc}/
for((i=1; i<=50; i++))
do
rm run$i/positions1.dump
rm run$i/positions2.dump
#rm run$i/nohup.out
done
echo "Removed position files for ${arc}!"
cd ..
#done
