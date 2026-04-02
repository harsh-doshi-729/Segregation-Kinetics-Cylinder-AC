# This script is used to plot the local monomer distribution of each region about the center of mass of the region
# This distribution is calculated by the script `local_monomer_distribution.c`

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import pathlib
from numpy.typing import ArrayLike
# importing the scripts that contains all the local paths
# The path to the directory containing this module is set as the environment variable PYTHONPATH
import system_file_paths as sysPaths 
import Segregation_Parameters as segParam
import PlottingTools as pt
# importing the regions config file:
sys.path.append(f"{sysPaths.POLYMER_PHYSICS}Scripts/Config_Files/")
import regions_config as reg
import plotMonomerDensity

# polymer information:
numberOfPolymers = 2
indexOfPolymer = 2
architecture = ""
numberOfRuns = 50
runIndex = -1

numberOfMonomers = 200 # each 

boxLength = 0 # The length of the cylindrical box; used to normalize the distance along the long axis while plotting
chosenDirection = "z" # default value for the direction letter

# Directory information:
acceptedDirectoryLabels = ["new_segregation", "Create_Initial_States"]
baseDirectory = sysPaths.NEW_SEGREGATION # The base directory where the data is stored

def SetBaseDirectory(directoryLabel: str) -> None:
    """Sets the base directory where the data is stored based on the directory label passed as an argument"""
    global baseDirectory
    if directoryLabel not in acceptedDirectoryLabels:
        print(f"The directory label '{directoryLabel}' is not accepted! Accepted labels are: {acceptedDirectoryLabels}")
        sys.exit(1)
    elif directoryLabel == acceptedDirectoryLabels[1]:
            baseDirectory = sysPaths.CREATE_INITIAL_STATES

def SetConstants():
    """Reads the argument(s) passed while invoking the script from the command line and sets it to the global variables"""
    global architecture
    global numberOfMonomers
    global runIndex
    global baseDirectory
    global boxLength
    global chosenDirection

    numberOfMandatoryArguments = 3
    if len(sys.argv) < numberOfMandatoryArguments + 1:
        print("Not enough arguments specified while invoking the script! Please pass the number of monomers, the name of the architecture, and the direction letter ('x', 'y', or 'z') in the following format:")
        print("python /<path>/plot.py <noOfMonomers> <architecture> <direction-letter>")
        print("Optionally, arguments for the run number and/or the directory label can be passed after these three arguments.")
        print("The run index is used to select the run from which the data is to be plotted. The directory label is used to select the directory from which the data is to be read.")
        print(f"Accepted directory labels are: {acceptedDirectoryLabels}")
        quit() # terminating the script
    
    architecture = sys.argv[2]
    # checking if the number of monomers and regions are valid:
    try:
        numberOfMonomers = int(sys.argv[1])
    except ValueError:
        print("The first argument could not be converted to an integer!")
        sys.exit(1)

    chosenDirection = sys.argv[3]
    if chosenDirection not in ['x', 'y', 'z']:
        print(f"The direction letter '{chosenDirection}' passed could not be recognized. Please pass one of 'x', 'y', and 'z'.")
        sys.exit(1)

    diameter = pt.ReadDiameter(numberOfMonomers, architecture)
    boxLength = segParam.ASPECT_RATIO * diameter

    if len(sys.argv) == numberOfMandatoryArguments + 1 + 1: # A singe optional argument is passed
        optionalArg = sys.argv[numberOfMandatoryArguments+1]
        if not optionalArg.isdigit(): # not a valid integer; must be the directory label
            SetBaseDirectory(optionalArg)
        else:
            runIndex = int(optionalArg)
    elif len(sys.argv) > numberOfMandatoryArguments + 1 + 1: # Two or more optional arguments are passed
        runIndex = int(sys.argv[numberOfMandatoryArguments+1])
        label = sys.argv[numberOfMandatoryArguments+2]
        SetBaseDirectory(label)

def GetDistributionFilePath(numberOfMonomers: int, architecture: str, runIndex: int, regionID: int) -> str:
    """
    Returns the path to the file containing the local monomer distribution data for the specified number of monomers, architecture, runIndex, and region ID.
    Args:
        numberOfMonomers (int): The number of monomers in the polymer.
        architecture (str): The architecture of the polymer.
        runIndex (int): The index of the run from which the data is to be read.
        regionID (int): The ID of the region for which the distribution is to be plotted, indexed from 1.
    """
    folderPath = sysPaths.GetFolder(baseDirectory, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)
    prefixFolder = ""
    if baseDirectory == sysPaths.CREATE_INITIAL_STATES:
        prefixFolder = "monomer_distribution/"
    return f"{folderPath}{architecture}/run{runIndex}/{prefixFolder}{chosenDirection}_local_monomer_distribution_reg{regionID}.csv"

def PlotLocalDistribution(ax: mpl.axes.Axes, numberOfMonomers: int, architecture: str, runIndex: int, regionID: int, legendLabel: str = "_ignored") -> None:
    """
    Reads the local monomer distribution data from the CSV file and plots it for a single region to the passed axis.
    Args:
        ax (mpl.axes.Axes): The axis on which to plot the local monomer distribution.
        numberOfMonomers (int): The number of monomers in the polymer.
        architecture (str): The architecture of the polymer.
        runIndex (int): The index of the run from which the data is to be read, indexed from 1.
        regionID (int): The ID of the region for which the distribution is to be plotted, indexed from 1.
    """
    df = pd.read_csv(GetDistributionFilePath(numberOfMonomers, architecture, runIndex, regionID))
    # The first column is the distance from the center of mass of the region
    # The second column is the local probability density of the monomers
    binLeftEdges = np.array(df.iloc[:, 0])
    binCentres = plotMonomerDensity.ShiftBinEdges(binLeftEdges)
    if segParam.RESCALE_LENGTHS:
        binCentres = binCentres / boxLength
    localDensity = np.array(df.iloc[:, 1])
    
    ax.plot(binCentres, localDensity, label = legendLabel, linestyle = "--", marker = ".")
    # Debugging: calculating area under the curves; should be equal to 1
    #areaUnderCurve = plotMonomerDensity.CalculateAreaUnderCurve(localDensity, binLeftEdges[1] - binLeftEdges[0])
    #print(f"Area under the curve for region {regionID} in run {runIndex} for architecture {architecture}: {areaUnderCurve:.4f}")

def ReadAndPlotDistribution(runIndex: int, showPlot: bool = False) -> None:
    """
    Reads the local monomer distribution data from the CSV file and plots it for each region.
    """
    if reg.USE_REGIONS:
        layout = (reg.NUMBER_OF_POLYMER_REGIONS, numberOfPolymers)
    else:
        layout = (1, numberOfPolymers)

    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/subplots_bold.mplstyle")
    fig, ax = plt.subplots(nrows = layout[0], ncols = layout[1], sharex = True, sharey = True)
    
    for regionID in range(reg.NUMBER_OF_REGIONS):
        ax_ij = pt.GetIndexedAxis(ax, regionID, layout)
        # Plotting:
        PlotLocalDistribution(ax_ij, numberOfMonomers, architecture, runIndex, regionID + 1)
        # ax_ij.legend(loc = 'upper right')
        ax_ij.set_title(f"  {reg.REGION_LABELS[regionID]}   ", loc = 'right', y = 0.9)
    
    xLabelModifier = ""
    # Rescaling the bin centres by box length
    if segParam.RESCALE_LENGTHS:
        xLabelModifier = r"/ $L$"
    fig.suptitle(f"Local Monomer Distribution along {chosenDirection} axis\n{sysPaths.SPECIAL_SIMULATION}, {architecture}, Run {runIndex}")
    fig.supxlabel(rf"Distance (${chosenDirection}$) from CoM of region {xLabelModifier}")
    fig.supylabel(rf"$P({chosenDirection}${xLabelModifier}$)$")

    if showPlot:
        plt.show(block = True)
    
    # Saving figure:
    saveFolder = f"{sysPaths.GetFolder(baseDirectory, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)}{architecture}/Analysis/Local_Monomer_Distribution/"
    path = pathlib.Path(saveFolder)
    path.mkdir(parents = True, exist_ok = True)  # Create the directory if it doesn't exist
    fig.savefig(f"{saveFolder}{chosenDirection}_local_monomer_distribution_run{runIndex}{segParam.FIG_EXT}")
    plt.close(fig)
    print(f"Local monomer distribution along {chosenDirection} direction plotted for {architecture} Run {runIndex}.")

def PlotComparison(runIndex: int, regionsOfInterest: list[int] = list(range(1, reg.NUMBER_OF_REGIONS + 1)), showPlot: bool = False) -> None:
    """
    Plots the comparison of local monomer distribution between multiple architectures for the given regions of interest.
    Args:
        runIndex (int): The index of the run from which the data is to be read, indexed from 1.
        regionsOfInterest (list[int]): A list of region IDs for which the comparison is to be plotted, indexed from 1.
    """
    layout = (1, len(regionsOfInterest))
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/subplots_bold.mplstyle")
    fig, ax = plt.subplots(nrows = layout[0], ncols = layout[1], sharex = True, sharey = True)
    counterIndex = 0

    for regionID in regionsOfInterest:
        ax_ij = pt.GetIndexedAxis(ax, counterIndex, layout)
        counterIndex += 1
        # Plotting:
        for architecture in segParam.AOI_COMPARE:
            PlotLocalDistribution(ax_ij, numberOfMonomers, architecture, runIndex, regionID, legendLabel = architecture)
        ax_ij.set_title(f"{reg.REGION_LABELS[regionID - 1]}   ", loc = 'center', y = 1)
        ax_ij.legend()
    xLabelModifier = ""
    # Rescaling the bin centres by box length
    if segParam.RESCALE_LENGTHS:
        xLabelModifier = r"/ $L$"
    fig.suptitle(f"Local Monomer Distribution along {chosenDirection} axis\n{sysPaths.SPECIAL_SIMULATION}, {architecture}, Run {runIndex}")
    fig.supxlabel(rf"Distance (${chosenDirection}$) from CoM of region {xLabelModifier}")
    fig.supylabel(rf"$P({chosenDirection}${xLabelModifier}$)$")

    if showPlot:
        plt.show(block = True)
    
    # Saving figure:
    saveFolder = f"{sysPaths.GetFolder(baseDirectory, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)}{architecture}/Analysis/Local_Monomer_Distribution/"
    path = pathlib.Path(saveFolder)
    path.mkdir(parents = True, exist_ok = True)  # Create the directory if it doesn't exist
    fig.savefig(f"{saveFolder}{chosenDirection}_local_comparison_run{runIndex}{segParam.FIG_EXT}")
    plt.close(fig)
    print(f"Comparison of local monomer distribution along {chosenDirection} axis plotted for {segParam.AOI_COMPARE} for Run {runIndex}.")


if __name__ == "__main__":
    SetConstants()
    # Testing SetConstants function:
    # print(f"Architecture: {architecture}, Number of Monomers: {numberOfMonomers}, Run Index: {runIndex}, Base Directory: {baseDirectory}")
    regionsOfInterest = [2, 4]

    if runIndex == -1: # Running for all runs
        for run in range(1, segParam.NUMBER_OF_RUNS + 1):
            if segParam.PLOT_COMPARISON:
                PlotComparison(run, regionsOfInterest)
            else:
                ReadAndPlotDistribution(run)
    else:
        if segParam.PLOT_COMPARISON:
            PlotComparison(runIndex, regionsOfInterest, showPlot = True)
        else:
            ReadAndPlotDistribution(runIndex, showPlot = True)




