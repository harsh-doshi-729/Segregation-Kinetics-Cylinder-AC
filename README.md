# Segregation Kinetics in Cylinder

This repository contains the codes and some representative data from the work in "Kinetics of segregation of topologically-modified ring polymers in cylindrical confinement".

There are four subfolders in this repository:
1. Create_Initial_States/ contains the codes and data for simulations that create the initial mixed state of the two polymers in cylindrical confinement.
2. Check_Concatenation/ contains the codes and data for simulations that check whether the mixed state created in Create_Initial_States/ has any unintentional concatenations between or within its polymers.
3. Segregation/ contains the codes and data for simulations where the polymers in the mixed state segregate.
4. Global_Scripts/ contains some common scripts that are further used in the other three subfolders.

Each of the three folders Create_Initial_States/, Check_Concatenation/, and Segregation/ has a subfolder Scripts/ that contains the codes (LAMMPS, C, Python) for that folder. Each of the three folders also has subfolders called "b200/" and "b500/". These folders contain data files about the polymers with 200 and 500 monomers respectively.
A sample simulation has been run for the Arc-1-2 architecture with the COM-Recenter initialization procedure. The data files for this simulation can be found in all three folders under the respective b200/fene_recenter/Arc-1-2/ subfolder.

Each of the three folders has a README.txt file that describes how the simulations were run in that folder.

## Using these codes:
This repository can be cloned and the codes can be directly run. 

Note that certain dependencies need to be met for the codes to run properly: LAMMPS, gcc (C compiler), and the Python libraries matplotlib, numpy, pandas, pathlib must be installed.
