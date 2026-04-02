# This script plots the segregation times distribution for an architecture but of differing special simulation cases. 
# By default, the comparison is with the special simulation and the fbsr_parallel case.
# Prerequisite: The special simulation variables in the config files must be unset
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import os
from typing import Dict, Tuple, List
import numpy as np
import numpy.typing as npt

#importing the scripts that contains all the local paths
# The path to the directory containing this module is set as the environment variable PYTHONPATH
import system_file_paths as sysPaths 
import Segregation_Parameters as segParam
import PlottingTools as pt

# importing the plotSegTimeDIstribution to resuse methods
sys.path.append(f"{sysPaths.NEW_SEGREGATION}Analysis/")
import plotSegTimeDistribution

# Global Variables:
numberOfMonomers = 200
architecture = "Arc0"
special_simulation_folders = [] # list to store all the special simulation folder
plottingStyle = "both"

numberOfSpecialSimulations = 1 # apart from the default one

def GetSpecialSimulationFolder(folderName: str) -> str:
    """Returns the absolute folder path corresponding to the special simulation folder name passed."""
    if folderName == 'default':
        return f"{sysPaths.NEW_SEGREGATION}b{numberOfMonomers}/"
    else:
        return sysPaths.GetFolder(sysPaths.NEW_SEGREGATION, numberOfMonomers, folderName)

def VerifySpecialSimulationFolder(special_simulation_folder: str) -> None:
    """Checks whether the simulation folder name passed exists in the NEW_SEGREGATION directory. Quits if it doesn't."""
    # First checking if the NEW_SEGREGATION directory is a default one:
    lastFolder = sysPaths.NEW_SEGREGATION.split(sep = '/')[-2] # extracting the folder before the last '/'
    if not (lastFolder == "new_segregation"):
        print("The file path for NEW_SEGREGATION is not the default one! Please check the System_File_Paths/system_file_paths.py config file.")
        sys.exit(1)
    
    expectedFolderPath = GetSpecialSimulationFolder(special_simulation_folder)
    if not os.path.isdir(expectedFolderPath):
        print(f"The directory for the special simulation folder {expectedFolderPath}/ could not be found!")
        sys.exit(1)

def VerifySpecialSimulationFolders(foldersList: List[str]) -> None:
    """Accepts a list of special simulations and verifies each one of them"""
    for folder in foldersList:
        VerifySpecialSimulationFolder(folder)

def AddDefaultFolder(folderList: List[str]) -> List[str]:
    """Returns a list with the added default folder if it does not exist already"""
    for folder in folderList:
        if folder == 'default':
            return folderList
    folderList.insert(0, 'default') # Insert element at the beginning of the list
    return folderList
    

def SetConstants():
    """Reads the command line arguments passed and sets the values of the global variable"""
    global numberOfMonomers
    global architecture
    global numberOfSpecialSimulations
    global special_simulation_folders
    global plottingStyle

    numberOfMandatoryArguments = 4 # atleast 3 arguments expected
    if len(sys.argv) < numberOfMandatoryArguments + 1:
        print("Not enough arguments passed! Please pass the number of monomers, architecture name, the number special simulations desired, the names of the special simulation folders as command line arguments.")
        print("For example: python /path/to/script/plotSegTimeComparison.py 200 Arc0 2 long_axis fbsr_parallel")
        print("Pass 'default' as a special simulation if you want to plot the default simulation data.")
        print(f"Optionally, an argument for the plotting style can be passed after these arguments. Available options: {list(plotSegTimeDistribution.plottingStyles.keys())}")
        sys.exit(1)

    architecture = sys.argv[2]
    # checking if the number of monomers is valid:
    try:
        numberOfMonomers = int(sys.argv[1])
    except ValueError:
        print(f"The entered numberOfMonomers {numberOfMonomers} cannot be converted to a number! Please provide a valid number")
        sys.exit(1)

    # Checking if the number of special simulations is valid:
    try:
        numberOfSpecialSimulations = int(sys.argv[3])
    except ValueError:
        print(f"The entered number of special simulations {numberOfSpecialSimulations} cannot be converted to a number! Please provide a valid number")
        sys.exit(1)
    if numberOfSpecialSimulations < 1:
        print("The number of special simulations must atleast be 1. Enter a number grater than equal to 1 and try again.")
        sys.exit(1)

    # reading special simulation folders:
    # The index of the first folder = numberOfMandatoryArguments
    # The index of the last folder = numberOfMandatoryArguments + numberOfSpecialSimulations - 1
    for i in range(numberOfMandatoryArguments, numberOfMandatoryArguments + numberOfSpecialSimulations):
        special_simulation_folders.append(sys.argv[i])
    # debugging:
    # print(f"Special Simulation folders: {special_simulation_folders}")
    VerifySpecialSimulationFolders(special_simulation_folders)
    
    # Adding the default folder:
    if numberOfSpecialSimulations == 1: # Only a single special folder passed
        special_simulation_folders = AddDefaultFolder(special_simulation_folders)
        if len(special_simulation_folders) <= 1: # checking if the list length increased
            print("Only the default folder name passed! Please enter more simulation folders to plot comparisons.")
            print("Alternatively, a single simulation folder apart from default can be passed. It will be compared with the default.")
            sys.exit(1)

    # Reading optional argument
    # For no optional arguments' case: total arguments = numberOfMandatoryArguments + numberOfSpecialSimulations -1 + 1
    if len(sys.argv) > numberOfMandatoryArguments + numberOfSpecialSimulations:
        plottingStyle = sys.argv[numberOfMandatoryArguments+numberOfSpecialSimulations]
        plotSegTimeDistribution.VerifyPlottingStyle(plottingStyle)

    # Updating number of special simulations:
    numberOfSpecialSimulations = len(special_simulation_folders)

def GetFolderPaths(foldersList: List[str]) -> Dict[str, str]:
    """Returns a dictionary of absolute folder paths with the key as the corresponding simulation folder name.
    The folder paths are computed for each of the folder names in the passed list"""
    folderPaths = {} # Empty dictionary?
    for folderName in foldersList:
        folderPaths[folderName] = GetSpecialSimulationFolder(folderName) # Adding to the dictionary
    return folderPaths

def ReadSegregationTimes(architecture: str) -> tuple[List[dict[int, int]], str]:
    """Reads the segregation times from the csv files of the comparison locations and returns them as a list
    Also returns the segregation criterion used for the comparison"""
    segTimesList = []
    folderPathDict = GetFolderPaths(special_simulation_folders)
    segregation_criterion = ""
    folderPathSuffix = f"{architecture}/Analysis/"

    for folderName in special_simulation_folders:
        segTimes, criterion = plotSegTimeDistribution.ReadSegTimes(f"{folderPathDict[folderName]}{folderPathSuffix}")
        segTimesList.append(segTimes)
        # Checking segregation criterion:
        if criterion == None:
            print(f"The segregation criterion could not be read from the file at {folderPathDict[folderName]}!")
            sys.exit(1)
        if len(segregation_criterion) == 0:
            segregation_criterion = criterion
        else:
            if segregation_criterion != criterion:
                print(f"All the segregation criteria for the segregation times are not the same! {segregation_criterion} and {criterion} cannot be compared")
                sys.exit(1)
    return segTimesList, segregation_criterion


def PlotComparison(showPlot: bool = False) -> None:
    """Plots and saves the comparison plots of segregation times for the default and special simulations."""
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/subplots_bold.mplstyle")

    segTimesList, segregation_criterion = ReadSegregationTimes(architecture)
    
    # Debugging: Printing special simulation folders:
    print(f"Special simulation folders: {special_simulation_folders}")
    # Titles for the plots:
    plotTitles = []
    for i in range(len(special_simulation_folders)):
        plotTitles.append(f"{special_simulation_folders[i]} ({len(segTimesList[i])} runs)")
    # Save file name:
    folder = f"{sysPaths.NEW_SEGREGATION}/Analysis/comparison/b{numberOfMonomers}/{architecture}_{numberOfSpecialSimulations}_comparisons_"
    # Statistics label:
    labels = []
    for i in range(len(segTimesList)):
        segTimesList[i] = plotSegTimeDistribution.RescaleTimesDict(segTimesList[i], segParam.TAU_0) # normalizing with respect to tau_0
        labels.append(plotSegTimeDistribution.GetStatisticsLabel(list(segTimesList[i].values())))
        # Printing statistics:
        print(f"{plotTitles[i]}:")
        plotSegTimeDistribution.PrintStatistics(list(segTimesList[i].values()), architecture)
    # Plotting:
    plottingFunctions = {"scatter": plotSegTimeDistribution.PlotScatterPlot, "distribution": plotSegTimeDistribution.PlotHistogram}
    
    # Convert each dictionary of segregation times to a list of values
    segTime2DList = [list(segTimesDict.values()) for segTimesDict in segTimesList] # Converting to a 2D List to convert to scientific notation
    segTime2DList, axisLabel = pt.ConvertMultipleArraysToScientificNotation(segTime2DList)
    # Convert back to dictionary after scientific notation conversion
    for i in range(len(segTimesList)):
        counter = 0
        for key in segTimesList[i].keys():
            segTimesList[i][key] = segTime2DList[i][counter]
            counter += 1
    xlabels = {"scatter": "Run Index", "distribution": rf"Segregation Time ({axisLabel}$\tau_0$)"}
    ylabels = {"scatter": rf"Segregation Time ({axisLabel}$\tau_0$)", "distribution": "Frequency"}
    
    for key in plotSegTimeDistribution.plottingStyles.keys():
        if plotSegTimeDistribution.plottingStyles[key]: # the 'both' key will be rejected by default
            fig, axes = plt.subplots(ncols = len(segTimesList), sharey = True, sharex = True)
            for i in range(len(segTimesList)):
                plottingFunctions[key](segTimesList[i], axes[i], labels[i], segregation_criterion, False)
                axes[i].set_title(plotTitles[i])
                axes[i].legend()
            fig.suptitle(f"Comparison of Segregation Times for {segParam.GetAliasArchitecture(architecture)}, {numberOfMonomers} monomers in each polymer")
            fig.supxlabel(xlabels[key])
            fig.supylabel(ylabels[key])

            if showPlot:
                plt.show(block = True)
            plt.close(fig)

            fig.savefig(f"{folder}{key}.png")

def PlotComparisonForAllArchitectures(architectures: List[str], showPlot: bool = False) -> None:
    """Plots the comparisons for all architectures passed in a single subplot figure"""
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/subplots_bold.mplstyle")

    segTimesDictList = []
    plotTitles = []
    segTimesList = []
    for architecture in architectures:
        list, segregation_criterion = ReadSegregationTimes(architecture)
        segTimesDictList += list
        segTimesList.append(np.array(list(segTimesDictList[-1].values())))
        # Titles for the plots:
        for i in range(len(special_simulation_folders)):
            plotTitles.append(f"{segParam.GetAliasArchitecture(architecture)} {special_simulation_folders[i]} ({len(segTimesDictList[i])} runs)")
    
    # Statistics label:
    labels = []
    for i in range(len(segTimesDictList)):
        segTimesList[i] = segTimesList[i] / segParam.TAU_0 # normalizing with respect to tau_0
        labels.append(plotSegTimeDistribution.GetStatisticsLabel(segTimesList[i]))
    # Reconstructing the dictionaries:
    for i in range(len(segTimesDictList)):
        segTimesDictList[i] = dict(zip(segTimesDictList[i].keys(), segTimesList[i]))
    # Plotting:
    plottingFunctions = {"scatter": plotSegTimeDistribution.PlotScatterPlot, "distribution": plotSegTimeDistribution.PlotHistogram}

    segTimesList, axisLabel = pt.ConvertMultipleArraysToScientificNotation(segTimesList)
    xlabels = {"scatter": "Run Index", "distribution": rf"Segregation Time ({axisLabel}$\tau_0$)"}
    ylabels = {"scatter": rf"Segregation Time ({axisLabel}$\tau_0$)", "distribution": "Frequency"}

    # Save file name:
    folder = f"{sysPaths.NEW_SEGREGATION}/Analysis/comparison/b{numberOfMonomers}/{len(architectures)}Arcs_{numberOfSpecialSimulations}comparisons_"

    for key in plotSegTimeDistribution.plottingStyles.keys():
        if plotSegTimeDistribution.plottingStyles[key]: # the 'both' key will be rejected by default
            fig, axes = plt.subplots(nrows = len(architectures), ncols = numberOfSpecialSimulations, sharex = True, sharey = True)
            for i in range(len(architectures)):
                for j in range(numberOfSpecialSimulations):
                    index = i * numberOfSpecialSimulations + j # rwo major notation for index in 2D array
                    plottingFunctions[key](segTimesDictList[index], axes[i][j], labels[index], segregation_criterion, False)
                    axes[i][j].set_title(plotTitles[index])
                    axes[i][j].legend()
            fig.suptitle(f"Comparison of Segregation Times for {len(architectures)} architectures, {numberOfMonomers} monomers in each polymer")
            fig.supxlabel(xlabels[key])
            fig.supylabel(ylabels[key])

            if showPlot:
                plt.show(block = True)
            plt.close(fig)

            fig.savefig(f"{folder}{key}.png")

# TODO: Perform T-test here as well

# Script:
if __name__ == "__main__":
    SetConstants()
    # Plotting styles flags:
    plotSegTimeDistribution.SetPlottingStyleFlagValues(plottingStyle)
    # Plotting comparison:
    if architecture == 'all':
        PlotComparisonForAllArchitectures(segParam.AOI_LIST, True)
    else:
        PlotComparison(True)
    
        


    

