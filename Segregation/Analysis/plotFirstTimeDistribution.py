# A script to plot the distribution of the first passage times of segregation; 
# The time when the CoM distance first crosses the first threshold f for the first time

# Give the option to specify a special simulation and an architecture for each subplot
# Or plot it for all AOIs for a given special simulation

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import sys

# importing the scripts that contains all the local paths
# The path to the directory containing this module is set as the environment variable PYTHONPATH
import system_file_paths as sysPaths 
import Segregation_Parameters as segParam
import PlottingTools as pt
# Importing other python files related to plotting segregation times:
import plotArcSegTimes
import plotSegTimeDistribution

# Global variables:
numberOfMonomers = 200
special_simulation = ""
directory = sysPaths.NEW_SEGREGATION

def SetConstants() -> None:
    """Accepts the arguments from the command line and sets the values to global variables"""
    global numberOfMonomers
    global special_simulation

    numberOfMandatoryArguments = 2
    if len(sys.argv) < numberOfMandatoryArguments + 1:
        print(f"Not enough arguments passed! Please pass the number of monomers and the special simulation name while invoking the script in the format specified below:")
        print("python plotFirstTimeDistribution.py <noOfMonomers> <special_simulation>")
        print(f"Accepted simulation names: {plotArcSegTimes.GetAcceptedSpecialSimulationNames(numberOfMonomers)}")
        quit()
    
    try:
        numberOfMonomers = int(sys.argv[1])
    except ValueError:
        print(f"The passed argument {sys.argv[1]} could not be converted to an integer. Please enter a valid number.")
        quit()
    
    special_simulation = sys.argv[2]
    # Verifying if the special_simulation is in the accepted names:
    acceptedNames = plotArcSegTimes.GetAcceptedSpecialSimulationNames(numberOfMonomers)
    if not special_simulation in acceptedNames:
        print(f"The simulation name {special_simulation} could not be recognized. Please enter one from the following list:")
        print(f"{acceptedNames}")
        sys.exit(1)
    
def GetFolder(special_simulation: str, architecture: str) -> str:
    """Returns the path to the analysis folder where the segregation times are stored for a particular architecture"""
    prefix = sysPaths.GetFolder(directory, numberOfMonomers, special_simulation)
    return f"{prefix}{architecture}/Analysis/"    

def ReadFirstTimes(folder: str, seg_criterion: str) -> dict[int, int]:
    """Reads the first passage times stored in a file directly under the passed folder and with the appropriate seg_criterion tag in the filename.
    Returns the read times as values in a dictionary with the run index being the corresponding key."""
    filePath = f"{folder}firstPassageTimes_{seg_criterion}.csv"
    return plotSegTimeDistribution.read_csv_as_dict(filePath)

def PlotFirstTimeBoxPlots(architectures: list[str], showPlot: bool = False) -> None:
    """Plots the box plots of the first passage times for the various architectures passed (as a list).
    A common special_simulation and segregation criterion (as mentioned in Segregation_Parameters.py) is considered for all architectures."""
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/bold.mplstyle")
    fig, ax = plt.subplots()
    firstTimesList = [] # A list of lists; one element list contains the first passage times for the runs of one architecture
    for architecture in architectures:
        folder = GetFolder(special_simulation, architecture)
        # Reading First Passage Times:
        firstTimesDict = ReadFirstTimes(folder, segParam.CRITERION_STRING)
        firstTimes = np.array(list(firstTimesDict.values())) / segParam.TAU_0 # Converting to units in terms of TAU_0
        firstTimesList.append(firstTimes)
    
    # Converting times to scientific notation:
    firstTimesList, axisLabel = pt.ConvertMultipleArraysToScientificNotation(firstTimesList)

    # Plotting box plots:
    customLabels = segParam.GetAliasArchitectureList(architectures)
    xlabel = "Architecture"
    if segParam.USE_CUSTOM_LABELS:
        customLabels = segParam.GetCustomLabels(architectures)
        xlabel = segParam.CUSTOM_AXIS_LABEL
    boxDict = ax.boxplot(x = firstTimesList, tick_labels = customLabels)
    # Axis labels and title:
    ax.set_title(f"Box Plots of First Passage Times for multiple architectures\n {numberOfMonomers} monomers, {special_simulation} simulation, {segParam.CRITERION_STRING}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(rf"First Passage Times ({axisLabel}$\tau_0$)")

    # Adding flier information:
    if segParam.INCLUDE_FLIERS_INFORMATION:
        plotArcSegTimes.AddFliersText(boxDict, ax)

    if showPlot:
        plt.show(block = True)
    
    # Saving figure:
    folder = sysPaths.GetFolder(directory, numberOfMonomers, special_simulation)
    fig.savefig(f"{folder}{special_simulation}_box_first_passage_times{segParam.FIG_EXT}")
    plt.close(fig)
    
if __name__ == "__main__":
    SetConstants()
    PlotFirstTimeBoxPlots(segParam.AOI_LIST, True)