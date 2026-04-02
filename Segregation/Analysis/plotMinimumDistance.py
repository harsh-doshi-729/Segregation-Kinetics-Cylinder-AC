# This script is used to plot the minimum distance between two polymers (overlapped or otherwise)
# when taken in an infinite cylindrcal confinement
# This minimum distance trajectory is also used to find the time of segregation for the polymers

from pathlib import Path
import matplotlib as mpl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import Segregation_Parameters as segParam

from numpy.typing import ArrayLike
# importing the scripts that contains all the local paths
# The path to the directory containing this module is set as the environment variable PYTHONPATH
import system_file_paths as sysPaths
import PlottingTools as pt
import plotCoMDistribution as com

# Importing the region Config file:
sys.path.append(f"{sysPaths.POLYMER_PHYSICS}Scripts/Config_Files/")
import regions_config as reg

# polymer information:
numberOfPolymers = 2
numberOfMonomers = 200 # each 
architecture = ""
numberOfRuns = 50
runIndex = -1
useRegions = False # Whether to use regions defined in the regions_config file
regionIDs = [1, 2] # The IDs (indexed from 1) for the regions between which the minimum distance is calculated

def SetConstants() -> None:
    """Reads the argument(s) while invoking the script and sets the values of the global variables."""
    global numberOfMonomers
    global architecture
    global runIndex
    global useRegions

    numberOfMandatoryArguments = 3

    if len(sys.argv) < numberOfMandatoryArguments + 1:
        print("Not enough arguments specified while invoking the script! Please pass the name of the number of monomers, the architecture, and a flag to indicate whether regions are used in the following format:")
        print("python /<path>/script.py <noOfMonomers> <architecture> <flag>")
        print("If the flag is true/1, the two region IDs for which the minimum distance is calculated must be passed as well.")
        print("An optional argument for the run index can be passed to plot the minimum distance for only that run.")
        sys.exit(1)
    else:
        try:
            numberOfMonomers = int(sys.argv[1])
        except ValueError:
            print(f"ERROR: The argument {sys.argv[1]} could not be interpreted as a number of monomers!")
            sys.exit(1)
        
        architecture = sys.argv[2]

        flag = sys.argv[3].lower()
        if flag in ["true", "1", "yes"]:
            useRegions = True
            numberOfMandatoryArguments += 2 # Two more arguments needed for the region IDs
            if len(sys.argv) < numberOfMandatoryArguments + 1:
                print("ERROR: The region IDs were not specified! Please pass the two region IDs (indexed from 1) for which the minimum distance is to be calculated.")
                sys.exit(1)
            else:
                try:
                    regionIDs[0] = int(sys.argv[4])
                    regionIDs[1] = int(sys.argv[5])
                except ValueError:
                    print(f"ERROR: The region IDs {sys.argv[4]} and {sys.argv[5]} could not be interpreted as integers!")
                    sys.exit(1)
                if regionIDs[0] < 1 or regionIDs[0] > reg.NUMBER_OF_REGIONS or regionIDs[1] < 1 or regionIDs[1] > reg.NUMBER_OF_REGIONS:
                    print(f"ERROR: The region IDs must be between 1 and {reg.NUMBER_OF_REGIONS} (inclusive) as per the regions_config file.")
                    sys.exit(1)

        if len(sys.argv) > numberOfMandatoryArguments + 1: # Optional argument
            try:
                runIndex = int(sys.argv[numberOfMandatoryArguments+1])
            except ValueError:
                print(f"ERROR: The optional argument {sys.argv[numberOfMandatoryArguments+1]} could not be interpreted as an integer for the run index.")
                sys.exit(1)

def GetFolder(numberOfMonomers: int, architecture: str) -> str:
    """Returns the folder where the files for the particular architecture are stored."""
    folderPrefix = sysPaths.GetFolder(sysPaths.NEW_SEGREGATION, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)
    return f"{folderPrefix}{architecture}/"

def ReadMinimumDistances(numberOfMonomers: int, architecture: str, runIndex: int) -> tuple[ArrayLike, ArrayLike]:
    """Reads the minimum distance calculated already from the run folder.
    Returns the timeSteps and the minimum distance as a tuple of lists."""
    folder = GetFolder(numberOfMonomers, architecture)
    if useRegions:
        fileName = f"min_distance_reg_{regionIDs[0]}_{regionIDs[1]}.dat"
    else:
        fileName = "min_distance.dat"
    filePath = f"{folder}run{runIndex}/{fileName}"

    # Reading with pandas:
    df = pd.read_csv(filePath)
    timeSteps = np.array(df.iloc[: , 0])
    minDistance = np.array(df.iloc[: , 1])

    return (timeSteps, minDistance)
                
def FindSegregationTime(timeSteps: ArrayLike, minDistances: ArrayLike) -> tuple[int, int]:
    """Finds the time of segregation based on the segregation criteria imposed on the minimum distance trajectory.
    Returns a tuple containing the first passage time and the segregation time."""
    # Note: The segregation criterion is set in Segregation_Parameters.py

    firstPassageTime = -1 # The time at which the polymers first get separated
    segTime = -1 # The time at which the polymers are said to be segregated

    firstThreshold = segParam.INF_F_THRESHOLD # The min distance must cross this threshold to be considered segregated
    secondThreshold = segParam.INF_S_THRESHOLD # The average of the subsequent min distance trajectory must be above this threshold
    intervalLength = int(segParam.INTERVAL_LENGTH * len(timeSteps)) # The window length over which the average should be performed

    # Iterating over each timeStep:
    for index, distance in zip(range(0, len(timeSteps)), minDistances):
        if distance > firstThreshold: # Checking for first condition
            # Setting first passage time if not already set:
            if firstPassageTime == -1:
                firstPassageTime = timeSteps[index]
            
            # Checking for second condition:
            meanDistance = np.mean(minDistances[index: index + intervalLength])
            if meanDistance > secondThreshold: # satisfied; polymers segregated
                segTime = timeSteps[index]
                break
            else:
                print(f"Segregation rejected at step {timeSteps[index]} since average was below second threshold.")
    return firstPassageTime, segTime

def PlotMinimumDistance(timeSteps: ArrayLike, minDistance: ArrayLike, segTime: int|None = None, showPlot: bool = False) -> None:
    """Plots the passed minimum distances against the time steps."""
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/bold.mplstyle")
    fig, ax = plt.subplots()

    if not useRegions:
        if segTime == None:
            firstPassageTime, segTime = FindSegregationTime(timeSteps, minDistance)
        # Scaling timesteps and segregation time
    didPolymersSegregate = True
    if segTime == -1 or useRegions:
        segTime = -1
        didPolymersSegregate = False
    timeSteps, segTime, axisLabel = com.ScaleSegregationTimes(timeSteps, segTime)

    # Plotting:
    ax.plot(timeSteps, minDistance)
    distLabel = "polymers"
    regionLabel = ""
    if useRegions:
        distLabel = "regions"
        regionLabel = f"Regions {reg.REGION_LABELS[regionIDs[0]-1]} and {reg.REGION_LABELS[regionIDs[1]-1]}\n"


    ax.set_title(f"Minimum Distance between two {distLabel}\n{regionLabel}{architecture}, {numberOfMonomers} monomers, Run {runIndex}")
    ax.set_xlabel(rf"Time ({axisLabel}$\tau_0$)")
    ax.set_ylabel("Minimum Distance")

    # Vertical line:
    if not useRegions:
        if didPolymersSegregate:
            com.PlotVerticalLine(ax, segTime)
        ax.legend()

    if showPlot:
        plt.show(block = True)

    # Saving figure:
    folder = GetFolder(numberOfMonomers, architecture)
    saveDirectory = f"{folder}Analysis/MinDistance/"
    Path(saveDirectory).mkdir(parents = True, exist_ok = True)
    fileRegionLabel = ""
    if useRegions:
        fileRegionLabel = f"_reg_{regionIDs[0]}_{regionIDs[1]}"
    saveFilePath = f"{saveDirectory}min_distance{fileRegionLabel}_r{runIndex}{segParam.FIG_EXT}"
    fig.savefig(saveFilePath)
    plt.close(fig)
    print(f"Minimum Distance plotted for {regionLabel}{architecture} Run {runIndex}.")

def WriteSegregationTimes(segTimes:list[int]) -> None:
    """Writes the list of segregation times to a file."""
    folder = GetFolder(numberOfMonomers, architecture)
    segTimesFilePath = f"{folder}Analysis/segregationTimes_{segParam.INF_CRITERION_STRING}.csv"
    segTimesFile = open(segTimesFilePath, "w")
    segTimesFile.write(f"Run Index, Segregation Time steps ({segParam.INF_CRITERION_STRING})\n") # Header

    for runIndex in range(numberOfRuns):
        segTime = segTimes[runIndex]
        if(segTime == -1):
            print(f"Run {runIndex + 1} did not show segregation\n")
        else:
            segTimesFile.write(f"{runIndex + 1}, {segTime}\n")
    segTimesFile.close()

def PlotMinimumDistancesForAll() -> None:
    """Plots the minimum distances for all runs."""
    global runIndex

    segTimes = [] # list to store all of the segregation times
    for runIndex in range(1, numberOfRuns+1):
        # Reading minimum distances:
        timeSteps, minDistance = ReadMinimumDistances(numberOfMonomers, architecture, runIndex)
        # Finding segregation time:
        segTime = None
        if not useRegions:
            firstPassageTime, segTime = FindSegregationTime(timeSteps, minDistance)
            segTimes.append(segTime)
        # Plotting minimum distances:
        PlotMinimumDistance(timeSteps, minDistance, segTime = segTime)
    if not useRegions:
        WriteSegregationTimes(segTimes)
    print(f"Plotted Minimum Distances for {numberOfRuns} run(s) of {architecture}.")
    if useRegions:
        print(f"The minimum distances between regions {reg.REGION_LABELS[regionIDs[0]-1]} and {reg.REGION_LABELS[regionIDs[1]-1]} were plotted.")
    
if __name__ == "__main__":
    SetConstants()
    if runIndex == -1: # No run index argument passed
        PlotMinimumDistancesForAll()
    else:
        timeSteps, minDistance = ReadMinimumDistances(numberOfMonomers, architecture, runIndex)
        PlotMinimumDistance(timeSteps, minDistance, showPlot = True)
