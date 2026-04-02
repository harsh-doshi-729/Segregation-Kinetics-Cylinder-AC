from platform import architecture
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path
import matplotlib as mpl

# importing the scripts that contains all the local paths
# The path to the directory containing this module is set as the environment variable PYTHONPATH
import system_file_paths as sysPaths 
import Segregation_Parameters as segParam
import PlottingTools as pt

# Global variables:
numberOfMonomers = 200
architecture = ""
runIndex = -1
special_simulation = ""
directory = sysPaths.NEW_SEGREGATION
directories = {"new_segregation": sysPaths.NEW_SEGREGATION, "Create_Initial_States": sysPaths.CREATE_INITIAL_STATES}

def SetConstants() -> None:
    """Accepts the arguments from the command line and sets the values to global variables"""
    global numberOfMonomers
    global architecture
    global runIndex
    global directory


    numberOfMandatoryArguments = 3
    if len(sys.argv) < numberOfMandatoryArguments + 1:
        print(f"Not enough arguments passed! Please pass the number of monomers and the special simulation name while invoking the script in the format specified below:")
        print("python plotThermoEnergy.py <noOfMonomers> <architecture> <run number>")
        print(f"Optionally, an argument for the directory can be passed after the above arguments. Accepted options: {list(directories.keys())}")
        quit()
    
    try:
        numberOfMonomers = int(sys.argv[1])
    except ValueError:
        print(f"The passed argument {sys.argv[1]} could not be converted to an integer. Please enter a valid number.")
        quit()
    
    architecture = sys.argv[2]

    try:
        runIndex = int(sys.argv[3])
    except ValueError:
        print(f"The passed argument {sys.argv[3]} could not be converted to an integer. Please enter a valid number.")
        quit()
    
    # Optional argument:
    if len(sys.argv) > numberOfMandatoryArguments + 1:
        directory_option = sys.argv[numberOfMandatoryArguments+1]
        # Verifying if the directory option is valid:
        if not directory_option in directories.keys():
            print(f"The directory option {directory_option} could not be recognized. Please enter one from the following list:")
            print(f"{list(directories.keys())}")
            sys.exit(1)
        else:
            directory = directories[directory_option]

def GetReadFilePath(numberOfMonomers: int, architecture: str, runIndex: int) -> str:
    """Returns the file path to the thermo output data file"""
    folder = sysPaths.GetFolder(directory, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)
    fileName = "thermo_output.dat"
    return f"{folder}{architecture}/run{runIndex}/{fileName}"

def PlotEnergies(filePath: str) -> None:
    """Reads the various energies written in the file located at filePath and plots the time series"""
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/big_bold.mplstyle")
    # checking existence of path:
    if not Path(filePath).exists():
        print(f"The thermo file path {filePath} does not exist!")
        quit()

    df = pd.read_csv(filePath, sep = ' ', skiprows = 1, header = None) # Skipping the first comment line in file
    time = df.iloc[:, 0]
    WCA_energy = df.iloc[:, 2]
    bond_energy = df.iloc[:, 3]
    potential_energy = df.iloc[:, 4] # not minimized at 0
    kinetic_energy = df.iloc[:, 5]
    total_energy = df.iloc[:, 6]
    potential_energy_offset = np.array(WCA_energy) + np.array(bond_energy) # Potential energy but offset such that minimum is at 0

    # Rescaling time:
    scaled_time = time / segParam.TAU_0
    scaled_time, axisLabel = pt.ConvertToScientificNotation(scaled_time)

    fig, ax = plt.subplots()
    # ax.plot(scaled_time, WCA_energy, label = "WCA")
    # ax.plot(scaled_time, bond_energy, label = "FENE")
    ax.plot(scaled_time, potential_energy_offset, label = "PE")
    ax.plot(scaled_time, kinetic_energy, label = "KE")
    ax.plot(scaled_time, total_energy, label = "Total")

    if segParam.SHOW_TITLE:
        ax.set_title(f"Time Series of System Energy\n{numberOfMonomers} monomers, {sysPaths.SPECIAL_SIMULATION} {architecture}, Run {runIndex}")
    else:
        ax.set_title(segParam.GetPaperFiguresTitle(numberOfMonomers, architecture))
    ax.set_xlabel(rf"Time ({axisLabel}$\tau_0$)")
    ax.set_ylabel(r"Energy ($k_B T$)")
    ax.legend()
    fig.subplots_adjust(bottom = 0.15)

    plt.show(block = True)

    # Saving plot:
    saveFolderPath = f"{sysPaths.GetFolder(directory, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)}{architecture}/Analysis/Energy/"
    Path(saveFolderPath).mkdir(exist_ok = True)
    saveFilePath = f"{saveFolderPath}energy_r{runIndex}{segParam.FIG_EXT}"
    fig.savefig(saveFilePath)
    print(f"Plotted the thermo energy for {architecture} Run {runIndex}.")
    
if __name__ == "__main__":
    SetConstants()
    filePath = GetReadFilePath(numberOfMonomers, architecture, runIndex)
    PlotEnergies(filePath)