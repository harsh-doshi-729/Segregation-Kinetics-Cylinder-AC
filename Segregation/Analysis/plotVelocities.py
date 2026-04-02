# This script plots the velocities with respect to some distance as a time series

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import Segregation_Parameters as segParam
from pathlib import Path

from typing import Tuple, List
from numpy.typing import ArrayLike
# importing the scripts that contains all the local paths
# The path to the directory containing this module is set as the environment variable PYTHONPATH
import system_file_paths as sysPaths
import PlottingTools as pt
# Importing the region Config file:
sys.path.append(f"{sysPaths.POLYMER_PHYSICS}Scripts/Config_Files/")
import regions_config as reg
# Importing the COM Distance plotting script:
import plotCoMDistribution as plotDist

# Global variables:
numberOfPolymers = 2
numberOfMonomers = 200
architecture = ""
runIndex = -1

velocityMode = "" # The mode used in the calculate_velocity.c script that indicates which distance was differentiated
acceptedModes = ["relative-com", "comparison"]
plotComparison = False # A flag to indicate whether the velocities should be plotted as part of a comparison with the COM distance

# Regions information:
useConfigRegions = False # A flag to indicate whether the regions in regions_config.py are used
regionIndices = [1, 2] # A list to store the region indices for relative velocity calculation

def SetConstants() -> None:
    """
    Sets the global variable values after reading the inputs from the command line.
    """
    global numberOfMonomers
    global architecture
    global runIndex
    global velocityMode
    global useConfigRegions
    global regionIndices
    global plotComparison

    numberOfMandatoryArguments = 4
    if len(sys.argv) < numberOfMandatoryArguments + 1:
        print("Not enough arguments passed! Please pass the number of monomers (integer), the architecture name (string), the run index (integer), and the velocity mode (string) as arguments to the command line.")
        print(f"Accepted velocity modes: {acceptedModes}")
        print("For example: python plot.py 200 Arc1_10 1 relative-com")
        print("If the velocities are to be plotted for all runs, pass 'all' as the run index.")
        print("A pair of optional arguments for the pair of region indices (eg.: 1 2) between which the relative velocity is calculated can be passed after the above arguments")
        sys.exit(1)
    else:
        try:
            numberOfMonomers = int(sys.argv[1])
        except ValueError:
            print("Invalid input for number of monomers. Please enter an integer.")
            sys.exit(1)

        architecture = sys.argv[2]
        try:
            runIndex = int(sys.argv[3])
        except ValueError:
            if runIndex == 'all':
                runIndex = -1
            else:
                print("Invalid input for run index. Please enter an integer.")
                sys.exit(1)

        velocityMode = sys.argv[4]
        if not velocityMode in acceptedModes:
            print(f"The passed velocity mode '{velocityMode}' was not recognized. Please pass one of {acceptedModes}.")
        elif velocityMode == "comparison":
            velocityMode = "relative-com"
            plotComparison = True

        # Optional arguments
        if len(sys.argv) > numberOfMandatoryArguments + 1:
            if len(sys.argv) < numberOfMandatoryArguments + 3:
                print("Invalid input for region indices. Please enter two integers.")
                sys.exit(1)
            else:
                try:
                    regionIndices = [int(i) for i in sys.argv[5:]]
                except ValueError:
                    print("Invalid input for region indices. Please enter integers.")
                    sys.exit(1)
                useConfigRegions = True

def GetFolder() -> str:
    """
    Returns the folder path for the architecture folder
    """
    return f"{sysPaths.GetFolder(sysPaths.NEW_SEGREGATION, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)}{architecture}/"

def GetVelocitiesFilePath(runIndex: int, useConfigRegions: bool = False, regionIndices: List[int] = [1, 2]) -> str:
    """
    Returns the file path for the velocities data file
    """
    folder = GetFolder()
    if useConfigRegions:
        regionLabel = f"_reg_{regionIndices[0]}_{regionIndices[1]}"
    else:
        regionLabel = ""
    return f"{folder}run{runIndex}/velocity_{velocityMode}{regionLabel}_run{runIndex}.csv"

def PlotVelocities(ax: mpl.axes.Axes, filePath: str, directionIndex: int = 2, label: str = "") -> tuple[str, str]:
    """
    Plots the velocities as a function of time on the passed axis
    Args:
        - ax : The axis on which to plot the velocities
        - filePath : The path to the velocities data file
        - directionIndex : The index of the direction to plot (0 for x, 1 for y, 2 for z)
        - label : The label for the plot
    Returns:
        A tuple of axis labels for x and y axes indicating the units of the scaled quantities
    """
    # Load the velocities data from the CSV file
    try:
        df = pd.read_csv(filePath)
    except FileNotFoundError:
        print(f"ERROR: File not found: {filePath}\n Please ensure the velocity mode is correct.")
        sys.exit(1)
    timeSteps = df.iloc[:, 0].to_numpy()
    # Scaling timesteps:
    timeSteps = timeSteps / segParam.TAU_0
    timeSteps, x_axisLabel = pt.ConvertToScientificNotation(timeSteps)

    velocities = df.iloc[: , 1 + directionIndex].to_numpy()
    # Scaling velocities:
    velocities *= segParam.TAU_0 # Multiplying by Tau_0 since time is in denominator; distance is in sigma units
    velocities, y_axisLabel = pt.ConvertToScientificNotation(velocities)

    # Plotting:
    ax.plot(timeSteps, velocities, label=label)

    return (rf" ({x_axisLabel}$\tau_0$)", rf" ({y_axisLabel}$\sigma/\tau_0$)")

def PlotVelocitiesAndSave(runIndex: int, showPlot: bool = False) -> None:
    """
    Reads the velocity data from the CSV file corresponding to the global variables, plots it as a timeseries, and saves the figure.
    Args:
        - runIndex : The index of the run to plot
        - showPlot : Whether to show the plot interactively
    """
    filePath = GetVelocitiesFilePath(runIndex, useConfigRegions, regionIndices)

    # Plotting:
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/bold.mplstyle")
    fig, ax = plt.subplots()
    x_axisLabel, y_axisLabel = PlotVelocities(ax, filePath)
    # Setting title and axis labels:
    if useConfigRegions:
        regionLabel = f"Regions: {reg.REGION_LABELS[regionIndices[0]-1]} & {reg.REGION_LABELS[regionIndices[1]-1]}\n"
    else:
        regionLabel = ""
    ax.set_title(f"{velocityMode} velocities for {numberOfPolymers} polymers\n{regionLabel}{numberOfMonomers} monomers each, Run {runIndex}")
    ax.set_xlabel(f"Time{x_axisLabel}")
    ax.set_ylabel(f"Velocity{y_axisLabel}")

    if showPlot:
        plt.show(block = True)

    # Saving figure
    saveFolder = f"{GetFolder()}Analysis/Velocities/"
    Path(saveFolder).mkdir(parents=True, exist_ok=True) # Creating the folder if it does not exist
    fileRegionLabel = ""
    if useConfigRegions:
        fileRegionLabel = f"_reg_{regionIndices[0]}_{regionIndices[1]}"
    fig.savefig(f"{saveFolder}velocity_{velocityMode}{fileRegionLabel}_r{runIndex}{segParam.FIG_EXT}")
    plt.close(fig)

    # Printing sucess message
    print(f"{velocityMode} velocity time series plotted for {architecture} Run {runIndex}\n{regionLabel}")

def PlotCOMDistance(ax: mpl.axes.Axes, runIndex: int) -> float:
    """
    Plots the z_COM Distance for a particular run using the plotCoMDistribution.py machinery
    Args:
        - ax : The axis on which the distance should be plotted
        - runIndex : The index of the run to plot
    Returns the segregation time when the polymers are said to have segregated
    """
    plotDist.numberOfPolymers = numberOfPolymers
    plotDist.numberOfMonomers = numberOfMonomers
    plotDist.architecture = architecture
    plotDist.runIndex = runIndex

    plotDist.SetConstants()

    # Plotting:
    timeSteps, z_com_list = plotDist.ReadData(runIndex)
    firstPassageTime, segTime = plotDist.FindImprovedSegregationTime(z_com_list)
    # Coverting timesteps to units of TAU_0
    timeSteps, segTime, axisLabel = plotDist.ScaleSegregationTimes(timeSteps, segTime)
    plotDist.PlotCoMDistance(ax, timeSteps, z_com_list, rf"{axisLabel}$\tau_0$")
    ax.set_title("") # Clearing title
    ax.set_xlabel("") # Clearing z axis label
    return segTime


def PlotComparisonWithDistance(runIndex: int, showPlot: bool = False) -> None:
    """
    Plots the COM Distance along with the relative COM velocities for the same run in two plots one below each other.
    Args:
        - runIndex: The index of the run to plot
        - showPlot: Whether to show the plot interactively
    """
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/subplots_bold.mplstyle")
    mpl.rcParams['axes.labelsize'] = 16
    mpl.rcParams['axes.labelweight'] = 'bold'
    mpl.rcParams['axes.formatter.limits'] = (-1, 4)
    fig, ax = plt.subplots(nrows = 2, ncols = 1, sharex=True)
    
    # Plotting distance:
    segTime = PlotCOMDistance(ax[0], runIndex)

    # Plotting velocity:
    filePath = GetVelocitiesFilePath(runIndex, False, [1, 2])
    x_axisLabel, y_axisLabel = PlotVelocities(ax[1], filePath)
    ax[1].set_ylabel(f"Velocity{y_axisLabel}")

    if segTime > 0: # Polymers segregated; plotting vertical line
        ax[0].axvline(segTime, color='r', linestyle='--', label = "Segregation Time")
        ax[1].axvline(segTime, color = 'r', linestyle = '--')
    ax[0].legend()

    # Titles and axis Labels:
    fig.suptitle(f"COM Distance and Velocity for {numberOfPolymers} {architecture} polymers\n{numberOfMonomers} monomers each, Run {runIndex}")
    fig.supxlabel(f"Time{x_axisLabel}")

    if showPlot:
        plt.show(block = True)

    # Saving figure:
    saveFolder = f"{GetFolder()}Analysis/Velocities/"
    Path(saveFolder).mkdir(parents=True, exist_ok=True) # Creating the folder if it does not exist
    fig.savefig(f"{saveFolder}velocity_{velocityMode}_comparison_r{runIndex}{segParam.FIG_EXT}")
    plt.close(fig)

    # Printing sucess message
    print(f"{velocityMode} velocity and Distance COM time series plotted for {architecture} Run {runIndex}.")

if __name__ == "__main__":
    SetConstants()
    if runIndex == -1:
        for runIndex in range(1, segParam.NUMBER_OF_RUNS + 1):
            PlotVelocitiesAndSave(runIndex)
    else:
        if plotComparison:
            PlotComparisonWithDistance(runIndex, showPlot = True)
        else:
            PlotVelocitiesAndSave(runIndex, showPlot = True)
