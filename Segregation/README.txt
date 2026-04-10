This directory contains files and scripts used to study the segregation of two polymers starting from a mixed state.
The mixed state must have been created by running the simulation in the Create_Initial_States/ directory.

Within this directory, there are the following subdirectories:
Scripts/: Contains the scripts used to run the segregation simulation and calculate some quantities for analysis.
Analysis/: Contains the scripts to analyse and plot the quantities calculated in the above directory.
b200/: Contains segregation data and plots for various initialization procedures and various architectures; with each polymer containing 200 monomers.
b500/: Contains segregation data and plots for various initialization procedures and various architectures; with each polymer containing 500 monomers.

In order to run the segregation simulations, the following steps should be performed:

i) Go to Scripts/. Run the execute.sh file using bash with the required arguments.
The usage will be displayed when the file is run without arguments. This will launch the segregation simulation in the background.

ii) In Scripts/, run the execute_postProcessing.sh script using bash with the required arguments.
The usage will be displayed when the file is run without arguments.
This will compile Scripts/general_monomer_density.c, radial_monomer_density.c, and calculateCoMDistribution.c.
These scripts will be run to calculate various monomer densities over the segregation run and to calculate the polymer COM trajectories.

iii) Go to Analysis/. Here, the files plotMonomerDensity.py can be run to plot the average monomer densities of the simulation.
Be sure to set the constants FIRST_THRESHOLD, SECOND_THRESHOLD, and CRITERION_STRING in Analysis/Segregation_Parameters.py.
This sets the segregation criterion for the next step.
The plotCoMDistribution.py script can be run to plot the CoM trajectories and determine the segregation time (based on the segregation criterion). 

iv) In Analysis/, the scripts plotSegTimeDistribution.py and plotArcSegTimes.py can be run to plot the distribution of 
segregation times over multiple runs and for multiple architectures.
The segregation times used in the paper have been included in their relevant data subdirectories as CSV files.
They can be found in Segregation/b<N>/<initialization_procedure>/<architecture>/Analysis/ directory.
Be sure to set the constants AOI_LIST in Analysis/Segregation_Parameters.py which lists the architectures of interest (AOI) that will be compared.

Note that there is a script Segregation_Parameters.py in Analysis/. This script defines many parameters that
are used in Python scripts thoughout this project. For the most part, the parameter values do not need to be changed from their default values.
But in some cases, like for steps (iii) and (iv), they should be changed according to the desired output. 
Some parameters can also be changed to modify the aesthetics of the plots generated.

