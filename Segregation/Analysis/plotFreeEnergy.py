# This script is useful for debugging. It plots the distribution of the reaction coordinate for the segregation run after segregation has occurred
# This distribution can then be used to find the free energy wrt a reaction coordinate

from matplotlib.axis import Axis
from matplotlib.figure import Figure
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from decimal import Decimal
import sys
from pathlib import Path
from numpy.typing import ArrayLike
# importing the scripts that contains all the local paths
# The path to the directory containing this module is set as the environment variable PYTHONPATH
import system_file_paths as sysPaths 
import Segregation_Parameters as segParam
import PlottingTools as pt
# importing the regions config file:
sys.path.append(f"{sysPaths.POLYMER_PHYSICS}Scripts/Config_Files/")
import regions_config as reg

# polymer information:
numberOfPolymers = 2
indexOfPolymer = 2
architecture = ""
numberOfRuns = 50
runIndex = -1

numberOfMonomers = 200 # each 
boxLength = 0 # Used to rescale the reaction coordinate while plotting free energy

reactCoordLabel = list(segParam.FREE_ENERGY_AXIS_LABELS.keys())[0] # The label/filename for the free energy distrbution based on the reaction coordinate chosen
regionIndices = []

operationModes = ["distribution", "free_energy", "manual"]
chosenMode = 3

def SetConstants():
    """Reads the argument(s) passed while invoking the script from the command line and sets it to the global variables"""
    global architecture
    global numberOfMonomers
    global runIndex
    global reactCoordLabel
    global boxLength
    global chosenMode

    numberOfMandatoryArguments = 6
    if len(sys.argv) < numberOfMandatoryArguments + 1:
        print("Not enough arguments specified while invoking the script! Please pass the number of monomers, the name of the architecture, an operation mode for the plot, the reaction coordinate label, and the pair of region IDs in the following format:")
        print("python /<path>/plotMonomerDensity.py <noOfMonomers> <architecture> <mode> <label> <regionID_1> <regionID_2>")
        print(f"Accepted operation modes: {operationModes}")
        print("Optionally, an argument for the run number can be passed after the above arguments.")
        quit() # terminating the script
    
    architecture = sys.argv[2]
    # checking if the number of monomers and regions are valid:
    try:
        numberOfMonomers = int(sys.argv[1])
    except ValueError:
        print("The first argument could not be converted to an integer!")
        sys.exit(1)
    
    # Reading diameter and setting box length
    diameter = pt.ReadDiameter(numberOfMonomers, architecture)
    boxLength = segParam.ASPECT_RATIO * diameter
    mode = sys.argv[3]
    try:
        chosenMode = operationModes.index(mode)
    except ValueError:
        chosenMode = -1
        print(f"ERROR: The passed operation mode '{mode}' was not recognized. Please choose on of the following accepted modes: {operationModes}")
        sys.exit(1)

    reactCoordLabel = sys.argv[4]

    # Reading region indices:
    try:
        region1 = int(sys.argv[5])
        region2 = int(sys.argv[6])
    except ValueError:
        print(f"ERROR: One or both of the arguments ('{sys.argv[5]}' and '{sys.argv[6]}') for the pair of the region IDs could not be converted to integers.")
        sys.exit(1)
    if region1 <= 0 or region2 <= 0:
        print(f"ERROR: One or both of the region IDs ({region1} and {region2}) are not positive.")
        sys.exit(1)
    regionIndices.append(region1)
    regionIndices.append(region2)


    # Constructing free energy label:
    if reg.USE_REGIONS:
        regionLabel = "reg"
    else:
        regionLabel = "pol"

    # Optional argument:
    if len(sys.argv) > numberOfMandatoryArguments + 1:
        try:
            runIndex = int(sys.argv[numberOfMandatoryArguments+1])
        except ValueError: # must be the free energy label if not an integer
            print(f"ERROR: The optional argument '{sys.argv[numberOfMandatoryArguments+1]}' could not be converted to an integer for the run index.")
            sys.exit(1)

def GetFolder(architecture: str) -> str:
    """Returns the folder where the segregation data is stored"""
    prefix = sysPaths.GetFolder(sysPaths.NEW_SEGREGATION, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)
    return f"{prefix}{architecture}/"


def GetFileNamePrefix(reactionCoordLabel: str, regionIndices: list[int]) -> str:
    """
    Returns the file name for the distribution or free energy file
    Args:
    - reactCoordLabel: The label indicating the reaction coordinate used for the distribution or free energy calculation
    - regionIndices: A pair of region IDs (indexed from 1) indicating the two regions for which the reaction coordinate is calculated
    """
    if reg.USE_REGIONS:
        regionLabel = "reg"
    else:
        regionLabel = "pol"
    return f"{reactionCoordLabel}_{regionLabel}_{regionIndices[0]}_{regionIndices[1]}"
       

def GetFileName(reactionCoordLabel: str, regionIndices: list[int]) -> str:
    """
    Returns the file name for the distribution or free energy file
    Args:
    - reactCoordLabel: The label indicating the reaction coordinate used for the distribution or free energy calculation
    - regionIndices: A pair of region IDs (indexed from 1) indicating the two regions for which the reaction coordinate is calculated
    """
    fileNamePrefix = f"{GetFileNamePrefix(reactionCoordLabel, regionIndices)}_{operationModes[chosenMode]}"
    fileSuffix = ""
    if runIndex == -1:
        fileSuffix = "_all"
    return f"{fileNamePrefix}{fileSuffix}.csv"

def ShiftBinEdges(binLeftEdges):
    """Accepts the list of left edges of the bins and shifts them to return a list of the bin centres"""
    binWidth = binLeftEdges[1] - binLeftEdges[0]
    binCentres = np.array(binLeftEdges) + binWidth/2
    return binCentres

def CalculateAreaUnderCurve(distribution: ArrayLike, binWidth: float) -> float:
    """Calculates and returns the area under the distribution curve based on the values and the bin width."""
    area = 0
    # Area = Sum(distribution[i] * binWidth)
    area = np.sum(distribution) * binWidth
    return area

def PlotDistributionFromFile(ax: mpl.axes.Axes, filePath: str, label: str = "") -> None:
    """
    Plots the probability distribution by reading the distribution from the passed file path.
    The distribution is plotted in the passed axis.
    """
    try:
        df = pd.read_csv(filePath)
    except FileNotFoundError:
        print(f"ERROR: The file '{filePath}' could not be found! Please check if the file exists or if the correct reaction coordinate is passed.")
        sys.exit(1)

    binsLeftEdges = df.iloc[:, 0]
    binWidth = float(binsLeftEdges[1]) - float(binsLeftEdges[0])
    distribution = df.iloc[:, 1]

    binCentres = ShiftBinEdges(binsLeftEdges)
    # Rescaling lengths with box length:
    if segParam.RESCALE_LENGTHS:
        binCentres = binCentres / boxLength
    
    # plotting histogram:
    ax.plot(binCentres, distribution, linestyle = "--", marker = '.', label = label)
    # plt.bar(binsLeftEdges, distribution, width = binWidth, align = 'edge', alpha = 0.5, label = f'Polymer {indexOfPolymer}')

    # Debugging: Calculating area under the curve
    print(f" Distribution Area: {CalculateAreaUnderCurve(distribution, binWidth)}")

def PlotDistributionOnAxis(ax: mpl.axes.Axes, reactCoord: str, runNumber: int, label: str = "", archName: str|None = None) -> None:
    """
    Plots the probability distribution after reading the file containing the distribution of the reaction coordinate located in the segregation simulation run folder.
    Plots the distribution on the passed axis.
    Args:
        - ax: The axis on which the distribution is to be plotted
        - reactCoord: The labeled for the reaction coordination used for calculating the distribution
        - runNumber: The run index of the simulation for which the distribution is being plotted
        - label: (optional) The label to be used in the legend for the plotted distribution
        - archName: (optional) The architecture name for which the distribution is calculated
    """

    # Case for when the distribution is computed over all runs; no run index passed
    runLabel = f"run{runNumber}"
    if runNumber == -1:
        runLabel = f"run1"
    # Setting file path:
    if archName == None:
        archName = architecture # Global architecture name
    filePath = f"{GetFolder(archName)}{runLabel}/{GetFileName(reactCoord, regionIndices)}"
    # Debugging:
    # print(f"Debugging: Reading distribution from file: {filePath}, archName = {archName}")
    # Debugging: Calculating area under the curve; printing header
    print(f"{label} {reactCoord}", end = "")
    PlotDistributionFromFile(ax, filePath, label = label)



def PlotDistribution(reactCoord: str, runNumber: int|None = None, showPlot: bool = False) -> None:
    """
    Plots the probability distribution of the reaction coordinate after reading the file located in the segregation simulation run folder.
    Sets the title and axis labels.
    Args:
        - reactCoord: The labeled for the reaction coordination used for calculating the distribution
        - runNumber: (optional) The run index of the simulation for which the distribution is being plotted
        - showPlot: A flag to indicate whether the plot should be shown in the interactive mode
    """
    # Setting run index if none is passed:
    if runNumber == None:
        runNumber = runIndex

    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/bold.mplstyle")
    fig, ax = plt.subplots()
    PlotDistributionOnAxis(ax, reactCoord, runNumber)
    
    # Setting axis labels and titles:
    if runNumber == -1:
        runLabel = f"{segParam.NUMBER_OF_RUNS} runs"
    else:
        runLabel = f"Run {runNumber}"
    axisLabel = segParam.FREE_ENERGY_AXIS_LABELS[reactCoord]
    labelModifier = ""
    if segParam.RESCALE_LENGTHS:
        labelModifier = r" / $L$"
    regionLabel = "" # The label for the names of the regions used
    if reg.USE_REGIONS:
        regionLabel = f"Regions: {reg.REGION_LABELS[regionIndices[0]-1]} and {reg.REGION_LABELS[regionIndices[1]-1]}\n"
    ax.set_title(f"Distribution of {axisLabel} for two {architecture} polymers\n{regionLabel}{numberOfMonomers} monomers each, Total: {numberOfPolymers} polymer(s), {runLabel}")
    ax.set_xlabel(f"{axisLabel}{labelModifier}")
    ax.set_ylabel("Probability Density")

    if showPlot:
        plt.show(block = True)

    # Saving figure:
    saveFolder = f"{GetFolder(architecture)}Analysis/Probability_Distribution/"
    Path(saveFolder).mkdir(parents = True, exist_ok = True)
    if runNumber == -1:
        fileRunLabel = "all"
    else:
        fileRunLabel = f"r{runNumber}"
    fileNamePrefix = GetFileNamePrefix(reactCoord, regionIndices)
    fig.savefig(f"{saveFolder}{fileNamePrefix}_{operationModes[chosenMode]}_{fileRunLabel}{segParam.FIG_EXT}")
    plt.close(fig)

    print(f"The probability distribution plotted for {fileNamePrefix} for {architecture} for {runLabel}.")

def PlotDistributionComparison(arcList: list[str], reactCoord: str, runNumber: int|None = None, showPlot: bool = False) -> None:
    """
    Plots the comparison of probability distributions of the reaction coordinate for different architectures in the same axis.
    Args:
        - arcList: A list of architecture names for which the distributions are to be plotted
        - reactCoord: The labeled for the reaction coordination used for calculating the distribution
        - runNumber: (optional) The run index of the simulation for which the distribution is being plotted
        - showPlot: A flag to indicate whether the plot should be shown in the interactive mode
    """
    # Setting run index if none is passed:
    if runNumber == None:
        runNumber = runIndex

    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/bold.mplstyle")
    fig, ax = plt.subplots()

    # Plotting all architectures:
    for arc in arcList:
        PlotDistributionOnAxis(ax, reactCoord, runNumber, label = arc, archName = arc)

    # Setting axis labels and titles:
    if runNumber == -1:
        runLabel = f"{segParam.NUMBER_OF_RUNS} runs"
    else:
        runLabel = f"Run {runNumber}"
    axisLabel = segParam.FREE_ENERGY_AXIS_LABELS[reactCoord]
    labelModifier = ""
    if segParam.RESCALE_LENGTHS:
        labelModifier = r" / $L$"
    regionLabel = "" # The label for the names of the regions used
    if reg.USE_REGIONS:
        regionLabel = f"Regions: {reg.REGION_LABELS[regionIndices[0]-1]} and {reg.REGION_LABELS[regionIndices[1]-1]}\n"
    ax.set_title(f"Distribution of {axisLabel} for two polymers\n{regionLabel}{numberOfMonomers} monomers each, Total: {numberOfPolymers} polymer(s), {runLabel}")
    ax.set_xlabel(f"{axisLabel}{labelModifier}")
    ax.set_ylabel("Probability Density")
    leg = ax.legend()
    # Increasing legend line thickness:
    for line in leg.get_lines():
        line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())

    if showPlot:
        plt.show(block = True)

    # Saving figure:
    saveFolder = f"{GetFolder(architecture)}Analysis/Probability_Distribution/" # Saved in the global architecture name folder
    Path(saveFolder).mkdir(parents = True, exist_ok = True)
    if runNumber == -1:
        fileRunLabel = "all"
    else:
        fileRunLabel = f"r{runNumber}"
    fileNamePrefix = GetFileNamePrefix(reactCoord, regionIndices)
    fig.savefig(f"{saveFolder}{fileNamePrefix}_{operationModes[chosenMode]}_comp_{fileRunLabel}{segParam.FIG_EXT}")
    plt.close(fig)

    print(f"The comparision of probability distributions has been plotted for {fileNamePrefix} for {arcList} for {runLabel}.")

def PlotManualComparison() -> None:
    """
    Plots the comparison of probability distributions for a particular reaction coordinate for any two manually chosen cases.
    """
    filePath1 = "/home/harsh/ComputationalPhysics/Polymer_Physics/LAMMPS_runs/cluster-data/new_segregation/b200/Previous_Attempts/inf_recenter/Arc1_10/run1/region_overlap_pol_1_2_distribution_all.csv"
    filePath2 = "/home/harsh/ComputationalPhysics/Polymer_Physics/LAMMPS_runs/cluster-data/new_segregation/b200/Previous_Attempts/inf_fbsr/Arc1_10/run1/region_overlap_pol_1_2_distribution_all.csv"
    
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/bold.mplstyle")
    fig, ax = plt.subplots()

    print("inf_recenter ", end = "")
    PlotDistributionFromFile(ax, filePath1, label = " INF Recenter (Antiparallel)")
    print("inf_fbsr ", end = "")
    PlotDistributionFromFile(ax, filePath2, label = "INF Glued-Monomers (Parallel)")
    ax.set_title(f"Distribution of overlap distance for two Arc1-10 polymers\n{numberOfMonomers} monomers each, 50 runs")
    ax.set_xlabel(rf"Overlap Distance $z$")
    ax.set_ylabel("Probability Density")
    ax.legend()

    plt.show(block = True)



def GetFreeEnergyFilePath(runIndex: int, reactionCoordLabel: str) -> str:
    """Returns the file path to the file containing the calculated free energy
    The argument reactionCoordLabel is a label that is suffixed in the file name descibing the reaction coordinate used"""
    folder = GetFolder(architecture)
    if runIndex == -1: # In the Analysis folder
        return f"{folder}Analysis/FreeEnergy/free_energy_{reactionCoordLabel}.csv"
    else: # In a run index folder:
        return f"{folder}run{runIndex}/free_energy_{reactionCoordLabel}.csv"

def ReadNumberOfSnapshots(filePath: str) -> int:
    """Reads the first line of the free energy file located at the filePath 
    and parses the number of snapshots used to calculate the free energy."""
    numberOfSnapshots = -1
    with open(filePath, "r") as file:
        firstLine = file.readline().strip('\n')
        try:
            numberString = firstLine.split(':')[-1].strip()
            numberOfSnapshots = int(numberString)
        except ValueError:
            print(f"The first line in the free energy file may not contain the header contianing number of snapshots!")
            sys.exit(1)
    return numberOfSnapshots


def GetRunLabel(runIndex: int) -> str:
    """Returns the run label printed in the title of the plot according to the passed run number"""
    if runIndex == -1:
        return  f"{numberOfRuns} Runs"
    else:
        return f"Run {runIndex}"
    
def GetFreeEnergyProfileFolder() -> str:
    """Returns the path to the folder where the plotted free energy profile is to be saved"""
    return f"{GetFolder(architecture)}Analysis/FreeEnergy/"

def GetFreeEnergyProfilePath(reactCoordLabel: str, runIndex: int) -> str:
    """Returns the path where the plotted free energy profile should be saved"""
    if runIndex == -1: # Path without a run index label
        return f"{GetFreeEnergyProfileFolder()}free_energy_{reactCoordLabel}.png"
    else: # Path with a run index label
        return f"{GetFreeEnergyProfileFolder}free_energy_{reactCoordLabel}r{runIndex}.png"

def PlotFreeEnergy(label: str = "", runIndex: int = runIndex, showPlot: bool = False) -> None:
    """Reads and plots the free energy with respect to a reaction coordinate. 
    The reaction coordinates are specified either as the label argument or in the FREE_ENERGY_LABELS dictionary in Segregation_Parameters.py"""

    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/bold.mplstyle")
    fig, ax = plt.subplots()

    # Iterating over different reaction coordinates:
    reactCoordLabels = {}
    numberOfSnapshots = 0 # The number of snapshots considered to calculate the free energy
    if len(label) > 0:
        reactCoordLabels[label] = label
    else:
        reactCoordLabels = segParam.FREE_ENERGY_LABELS
    for label_key in reactCoordLabels.keys():
        freeEnergyFilePath = GetFreeEnergyFilePath(runIndex, label_key)
        # Reading the number of snapshots:
        numberOfSnapshots = ReadNumberOfSnapshots(freeEnergyFilePath)
        # Reading the free energy; after the first line
        df = pd.read_csv(freeEnergyFilePath, skiprows = 1) 
        reactionCoord = df.iloc[:, 0]
        stringFreeEnergy = df.iloc[:, 1] # Can it read the infinity values too?
        freeEnergy = [] # List of float values of the free energy
        # conveting to float from string:
        for string in stringFreeEnergy:
            freeEnergy.append(float(string))

        # Rescaling reaction coordinate with box length
        if segParam.RESCALE_LENGTHS:
            reactionCoord = np.array(reactionCoord) / boxLength

        # Plotting
        ax.plot(reactionCoord, freeEnergy, linestyle = "--", marker = '.', label = reactCoordLabels[label_key])
    xLabelModifier = ""
    if segParam.RESCALE_LENGTHS:
        xLabelModifier = r"/ $L$"
    runLabel = GetRunLabel(runIndex)
    ax.set_title(f"Free energy profile for system of {numberOfPolymers} {architecture} polymer(s)\n {numberOfMonomers} monomers each, {runLabel}, {numberOfSnapshots} snapshots")
    ax.set_xlabel("Distance between region CoMs along z axis %s" % (xLabelModifier))
    ax.set_ylabel("Free Energy (kT)")
    leg = ax.legend()
    # Increasing legend line thickness:
    for line in leg.get_lines():
        line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())

    if showPlot:
        plt.show(block = True)
    # saving figure:
    if len(label) > 0:
        label = label + '_'
    Path(GetFreeEnergyProfileFolder()).mkdir(parents = True, exist_ok = True)
    fig.savefig(GetFreeEnergyProfilePath(label, runIndex))
    plt.close(fig)

    print(f"Free Energy profile plotted for {runLabel} of {architecture}.")

if __name__ == "__main__":
    SetConstants()
    print(f"Chosen Mode: index = {chosenMode}, mode = {operationModes[chosenMode]}")
    if operationModes[chosenMode] == "distribution":
        if segParam.PLOT_COMPARISON:
            PlotDistributionComparison(segParam.AOI_COMPARE, reactCoordLabel, runIndex, True)
        else:
            PlotDistribution(reactCoordLabel, showPlot = True)
    elif operationModes[chosenMode] == "free_energy":
        PlotFreeEnergy(reactCoordLabel, runIndex, True)
    elif operationModes[chosenMode] == "manual":
        PlotManualComparison()
