This directory contains files and scripts used to create the initial states or configuration files from where segregation simulations can begin.
The initial states are usually mixed states of two polymers confined in a cylinder. They are created using one of the initialization procedure algorithms.

The b200/ and b500/ directories contain the files for 200 monomer and 500 monomer systems respectively.
These directories contain a separate directory for each of the initialization procedure used: fene_recenter/ (COM-Recenter), fene_glued/ (Glued-Monomers), fene_ladder/ (Ladder-Like), and fene_mutual/ (Mutual-Attraction).
The b200/ directory contains an additional directory inf_recenter/ which is a simulation of polymers in a nearly infinite cylinder.

In order to create a mixed state for a particular architecture, the following procedure is followed:
i) An initial state file is created for a particular polymer architecture which will be used to start the simulation.
ii) The input file is then fed into a LAMMPS script such as shrink_relax.lammps to run a simulation. The state obtained at the end of the simulation is dumped and is the required mixed state.

After the mixed states are created, they are tested to see if they are suitable to proceed to the segregation step. Currently these tests include:
i) Check for concatenations across the two polymers and within each polymer
ii) Check for density anomalies caused due to the shrink-relax algorithm
iii) Check if the two polymers are sufficiently mixed

Check (i) is done in the Check_Concatenation directory located at ../
Check (ii) is done by calculating monomer distributions in the cylinder after the shrinking-relax is complete. The relevant files are stored in the Create_Initial_States/b<N>/<initialization_procedure>/<architecture>/run1/ folder.
Check (iii) is carried out by performing a long simulation of the mixed state and observing the CoM distance. The relevant files are stored in
the Create_Initial_States/b<N>/<initialization_procedure>/<architecture>/run1/ folder.

These steps can be done at once sequentially one after the other using the bash script Create_Initial_States/Scripts/execute.sh.
The usage of this script is:
$ bash execute.sh <number-of-monomers> <architecture> <run-number> <seed> <initialization-procedure-folder>

Alternatively, the individual LAMMPS scripts are present in Create_Initial_States/Scripts/ and can be run manually by the user. These scripts take the following variables while being launched:
seed, numberOfMonomers, run (run-number), angle1 (orientation angle of 1st polymer), angle2 (orientation angle of the second polymer), and final_radius (the desired final radius of the cylinder after the shrinking procedure)
Example usage:
<LAMMPS_EXEC> -in fene_recenter.lammps -var seed 532905 -var numberOfMonomers 200 -var run 01 -var angle1 0 -var angle2 0 -var finalRadius 3.12 -log shrink_log.lammps
An initial_configuration.txt file like in Create_Initial_States/Scripts/ must be present in the directory where this script is run.
This script will perform the simulation and write the snapshots of the system as a LAMMPS output file. It also stores the last snapshot in a special langevin_mixed_state_01.txt file which is used for the segregation simulation.

The checks (ii) and (iii) are also performed automatically in execute.sh after the simulation has finished.
This analysis involves running the C files present in Create_Initial_States/Scripts/. They can be compiled and run manually as well. Please ensure that the path to the base directory (Segregation-Kinetics-Cylinder-AC/) is passed as an include path using the -I option while compiling any script using gcc.
The C scripts take some command line arguments, and the usage is explained when the compiled output is run without passing any arguments.

Finally, the Python scripts in Create_Initial_States/Scripts/ can be run manually to plot the calculated monomer densities that were also shown in the paper.
The scripts are prefixed with the word 'plot'. All these scripts take a few command line arguments to run properly. The usage is displayed when these scripts are run without passing any arguments.

Once all the above steps are done, the user can move on to the Check_Concatenation/ folder if the concantenations have not already been checked.

