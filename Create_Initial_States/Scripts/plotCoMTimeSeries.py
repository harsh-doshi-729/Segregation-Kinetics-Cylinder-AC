# This script plots the CoM of different regions in the system as a function of time during the shrink-relax equilabration phase

from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import ArrayLike
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
from pathlib import Path


# importing the scripts that contains all the local paths
sys.path.append(f"../../Global_Scripts/System_File_Paths/") # Adding the path to the system file paths module to the system path
import system_file_paths as sysPaths 
sys.path.append(f"{sysPaths.SEGREGATION}Analysis/")
import Segregation_Parameters as segParam
# Importing plotting tools:
sys.path.append(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/")
import PlottingTools as pt
# import plotCoMDistribution:
import plotCoMDistribution as plotCoM
# import regions_config file
sys.path.append(f"{sysPaths.GLOBAL_SCRIPTS}Config_Files/")
import regions_config as reg

# polymer information:
numberOfPolymers = 2
indexOfPolymer = 1
architecture = "Arc2"
runIndex = 1
boxLength = 20

numberOfRuns = 50
numberOfMonomers = 200 # each 
numberOfSteps = 10 ** 7
dataInterval = 100
numberOfBins = 65
numberOfMandatoryArguments = 3
readFileNamePrefix = "com"
initializationProcedure = "" # The name of the initialization procedure used to generate the initial state; this is used to read from the correct folder

TAU_0 = 200 # The time scale in no of iterations for the LJ units; inverse of the timestep chosen
# TODO: move this to SegregationParameters.py

def SetConstants():
    """Sets the global variables according to the case pertaining to the number of monomers used"""
    global boxLength
    global numberOfSteps
    global dataInterval

    diameter = pt.ReadDiameter(numberOfMonomers, architecture, initializationProcedure)
    boxLength = segParam.ASPECT_RATIO * diameter
    if numberOfMonomers == 200:
        numberOfSteps = 2 * 10 ** 7 # 5 times relaxation time of a simple ring polymer
        dataInterval = 1000
    elif numberOfMonomers == 500:
        numberOfSteps = 1.25 * 10 ** 8 # 5 times relaxation time of a simple ring polymer
        dataInterval = 1000

def SetArchitectureAndMonomers():
    """Reads the argument(s) passed while invoking the script from the command line and sets it to the architecture and numberOfMonomers"""
    global architecture
    global numberOfMonomers
    global initializationProcedure
    
    if len(sys.argv) < numberOfMandatoryArguments + 1:
        print("Not enough arguments specified while invoking the script! Please pass the number of monomers (integer), the architecture (string), and the initialization procedure (string) in the following format:")
        print("python /<path>/plotCoMDistribution.py <noOfMonomers> <architecture> <initializationProcedure>")
        print("An optional argument for the run index can be passed to plot the CoM data for only that run")
        print(f"Please ensure the Tau_0 value has been correctly set! Current value: {TAU_0}")
        quit() # terminating the script
    
    try:
        numberOfMonomers = int(sys.argv[1])
    # checking if the number of monomers is valid:
    except ValueError:
        print(f"The entered numberOfMonomers {numberOfMonomers} cannot be converted to a number! Please provide a valid number")
        sys.exit(1)

    architecture = sys.argv[2]
    initializationProcedure = sys.argv[3]

def GetFolder(baseFolder: str, numberOfMonomers: int, architecture: str) -> str:
    """Returns the architecture folder after taking into account the special simulation, if any."""
    return f"{sysPaths.GetFolder(baseFolder, numberOfMonomers, initializationProcedure)}{architecture}/"

def GetRunIndex():
    """Reads the third argument passed while invoking the script from the command line and returns it if it is an integer. This argument is optional"""
    if len(sys.argv) == numberOfMandatoryArguments + 1: # just file name, the architecture, and number of monomers
        return -1
    else:
        runIndex = sys.argv[numberOfMandatoryArguments + 1]
        if runIndex.isdigit():
            return int(runIndex)
        else:
            print(f"The third argument was expected to be the run index, but the one provided {runIndex} cannot be cast to an integer")
            exit(1)

def PlotTimeSeries(z_com, ax, indexOfPolymer: int, timeSteps = None):
    # plotting time series:
    if timeSteps is None:
        timeSteps = range(0, numberOfSteps, dataInterval)
    ax.plot(timeSteps, z_com, label = f"Polymer {indexOfPolymer}")
    ax.set_xlabel("Time Steps")
    ax.set_ylabel("z coordinate of CoM")
    ax.set_title(f"Time Evolution of centre of mass\n {numberOfPolymers} {architecture} polymers, Run {runIndex}")

def PlotDistribution(z_com, ax, indexOfPolymer: int):
    # plotting:
    ax.hist(z_com, numberOfBins, alpha = 0.5, label = f"Polymer {indexOfPolymer}")
    ax.set_title(f"Histogram of Z coord of CoM\n {numberOfPolymers} {architecture} polymers, Run {runIndex}")
    ax.set_xlabel("z coordinate")
    ax.set_ylabel("Frequency")

def PlotVerticalLine(axes, time: int):
    axes.axvline(x = time, color = 'r', linestyle = '--', label = 'Time of Segregation')


def ReadData(runIndex: int, readFileNamePrefix: str) -> tuple[ArrayLike, list[ArrayLike]]:
    """Reads the z coordinate of the CoM Data for both the polymers and returns the (timeSteps, double dimesnional list containing the read data) for a particular run.
        The readFileNamePrefix variable is the common prefix of the CoM data files. For example for the data files labelled as 'langevin_single_com1.dat', the prefix is 'langevin_single_com'"""
    folder = GetFolder(sysPaths.CREATE_INITIAL_STATES, numberOfMonomers, architecture)
    z_com_list = []
    timeSteps = []
    isTimeStepsSet = False # a flag to indicate whether the time steps has been read and set already
    limitingSteps = numberOfSteps # the minium number of steps in the two CoM files; using this to make the shapes of both z_com arrays the same

    for indexOfPolymer in range(1, numberOfPolymers+1):
        filePath = f"{folder}run{runIndex}/{readFileNamePrefix}{indexOfPolymer}.dat"
        # print(f"Reading the file: {filePath}") # Debugging
        df = pd.read_csv(filePath)
        z_com = np.array(df.iloc[: , 3])

        limitingSteps = min(limitingSteps, len(z_com))
        if not isTimeStepsSet:
            timeSteps = np.array(df.iloc[:, 0])
            isTimeStepsSet = True

        z_com_list.append(z_com)

    # making the shape of both lists the same; cutting out overhangs
    timeSteps = timeSteps[:limitingSteps]
    for index in range(numberOfPolymers):
        z_com_list[index] = z_com_list[index][:limitingSteps]

    return timeSteps, z_com_list

def PlotData(timeSteps, z_com_list, showPlot: bool = False):
    """Takes the timesteps and the z_com data of the two polymers and plots one against the other. Plots a vertical red line at the time = segTime"""
    folder = GetFolder(sysPaths.SEGREGATION, numberOfMonomers, architecture)
    # fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()

    # average:
    # z_com = (z_com1 + z_com2)/2
    # PlotDistribution(z_com, ax1, indexOfPolymer)
    for index in range(numberOfPolymers):
        PlotTimeSeries(z_com_list[index], ax2, index+1, timeSteps)

    # PlotCoMDistance(z_com_list[0], z_com_list[1], ax2, timeSteps)
    # ax1.legend()
    leg = ax2.legend()
    # Increasing legend line thickness:
    for line in leg.get_lines():
        line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())

    if showPlot:
        plt.show(block = True)
    # fig1.savefig(f"{folder}Analysis/CoM_Distribution/CoMDistribution_r{runIndex}.png")
    # fig2.savefig(f"{folder}Analysis/TimeSeries/CoMTimeSeries_r{runIndex}.png")
    
    #closing figures to save memory:
    # plt.close(fig1)
    plt.close(fig2)
    print(f"Plotted CoM Time Series for {architecture} Run {runIndex}")

def ReadAndPlotPolymerCoM(runIndex: int, showPlot: bool = False) -> None:
    """Reads the CoMs of individual polymers and plots the z CoM distance between them"""
    timeSteps, z_CoM_List = ReadData(runIndex, "com_reg") # Reading the com#.dat files
    # Scaling timesteps:
    scaledTimeSteps = timeSteps / TAU_0
    scaledTimeSteps, axisLabel = pt.ConvertToScientificNotation(scaledTimeSteps)

    folder = GetFolder(sysPaths.CREATE_INITIAL_STATES, numberOfMonomers, architecture)
    CoMFolder = f"{folder}Analysis/Polymer_CoM/"
    Path(CoMFolder).mkdir(parents = True, exist_ok = True)
    saveFilePath = f"{CoMFolder}CoM_Distance_r{runIndex}{segParam.FIG_EXT}"
    PlotCoMDistance(saveFilePath, runIndex, z_CoM_List, scaledTimeSteps, axisLabel, showPlot)
    # Plotting CoM Time series too:
    # PlotPolymerCoMsSeparately(f"{CoMFolder}CoM_TimeSeries_r{runIndex}.png", runIndex, z_CoM_List, scaledTimeSteps, axisLabel, showPlot)

def PlotPolymerCoMForAllRuns(showPlots: bool = False) -> None:
    """Plots the polymer CoM distance for all runs by directly reading polymer CoM files"""
    for runIndex in range(1, numberOfRuns + 1):
        try:
            # timeSteps, z_com_list = ReadData(runIndex, readFileNamePrefix)
            # PlotData(timeSteps, z_com_list, showPlots)
            ReadAndPlotPolymerCoM(runIndex, showPlots)
        except FileNotFoundError as fileError:
            print(f"The file {fileError.filename} could not be opened for run {runIndex}!")
            print(f"ERROR: {fileError.strerror}")
            print(f"Continuing to the rest of the script: {sys.argv[0]}")

def FindOrientation(z_CoM_List: list[ArrayLike]) -> str:
    """Find the relative orientation of the two polymers in the last snapshot of the run.
    Orientation can only be found for polymers containing 4 regions.
    Accepts the list of numpy arrays containing time series of CoMs for each region.
    Returns 'Parallel' if the polymers are in a parallel orientiation.
    Returns 'Antiparallel' if the polymers are in an antiparallel orientation.
    Returns 'Ambiguous' for all other cases."""

    # Calculating 1D vectors joining big loop to smaller loops:
    orientVectors = []
    for n in range(numberOfPolymers):
        bigRegID = n * reg.NUMBER_OF_POLYMER_REGIONS # The region ID for the big loop
        smallRegID = n * reg.NUMBER_OF_POLYMER_REGIONS + 1 # The region ID for the small looped region
        orientVectors.append(z_CoM_List[bigRegID][-1] - z_CoM_List[smallRegID][-1]) # calculating vector at the last step

    dotProduct = orientVectors[0] * orientVectors[1] # Only the z-component dot product
    # Classifying orientation:
    if dotProduct > 0: # positive = parallel
        return 'Parallel'
    elif dotProduct < 0: # negative = antiparallel
        return 'Antiparallel'
    else: # perpendicular? should never occur
        print(f"The relative orientation seems ambiguous for {architecture} Run {runIndex}. Dot Product = 0.")
        return 'Ambiguous'

def ReadAndPlotRegionCoM(runIndex: int, showPlot: bool = False) -> str:
    """Reads the CoM Time series data and plots them.
    Returns the oreintation (string) of the polymers: Parallel, Antiparallel, or Ambiguous"""
    folder = GetFolder(sysPaths.CREATE_INITIAL_STATES, numberOfMonomers, architecture)
    timeSteps, z_CoM_List = plotCoM.ReadRegionData(folder, runIndex)
    # scaling timesteps by the time it takes to generate 1 mixed state
    # so that the time is reported in terms of which run mixed state was dumped last
    # scaledTimeSteps = timeSteps / numberOfSteps

    # Scaling timesteps by tau:
    scaledTimeSteps = timeSteps / TAU_0 # type: ignore
    scaledTimeSteps, axisLabel = pt.ConvertToScientificNotation(scaledTimeSteps)
    # Scaling the CoM positions by box length
    yLabelModifier = ""
    if segParam.RESCALE_LENGTHS:
        z_CoM_List = pt.RescaleArraysByLength(z_CoM_List, boxLength)
        z_CoM_List, yAxisLabel = pt.ConvertMultipleArraysToScientificNotation(z_CoM_List)
        if len(yAxisLabel) != 0:
            yAxisLabel = f"({yAxisLabel})"
        yLabelModifier = rf"/ $L$ {yAxisLabel}"

    fig, ax = plt.subplots()

    for i in range(len(z_CoM_List)):
        ax.plot(scaledTimeSteps, z_CoM_List[i], label = reg.REGION_LABELS[i])
    if segParam.SHOW_TITLE:
        ax.set_title(f"CoM Time Series of various regions\n {numberOfPolymers} {architecture} polymer(s), {numberOfMonomers} monomers each, Run {runIndex}")
    else: # Printing a shorter title
        ax.set_title(segParam.GetPaperFiguresTitle(numberOfMonomers, architecture))
    ax.set_xlabel(rf"Time ({axisLabel}$\tau_0$)")
    ax.set_ylabel(r"$z_{CoM}$ %s" % (yLabelModifier))
    leg = ax.legend()
    # Increasing legend line thickness:
    for line in leg.get_lines():
        line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)
        
    if showPlot:
        plt.show(block = True)
    
    # Calculating polymer orientation:
    orientation = FindOrientation(z_CoM_List)
    print(f"{architecture} Run {runIndex}: {orientation} orientation")


    # saving figure:
    analysisFolder = f"{folder}Analysis/TimeSeries/"
    Path(analysisFolder).mkdir(parents = True, exist_ok = True)
    fig.subplots_adjust(bottom = 0.2, left = 0.2)
    fig.savefig(f"{analysisFolder}CoM_TimeSeries_r{runIndex}{segParam.FIG_EXT}")
    plt.close(fig)
    print(f"Plotted CoM Time Series for {architecture} Run {runIndex}")
    return orientation

def PlotRegionTimeSeriesForAllRuns(showPlots: bool = False):
    """Reads the region CoM time series and plots it for all runs
    Defaults to polymers as regions if no regions are used"""
    orientations = {"Parallel": [], "Antiparallel": [], "Ambiguous": []} # to store the list of runs for each orientation
    for runIndex in range(1, numberOfRuns + 1):
        try:
            # timeSteps, z_com_list = ReadData(runIndex, readFileNamePrefix)
            # PlotData(timeSteps, z_com_list, showPlots)
            orientation = ReadAndPlotRegionCoM(runIndex)
            orientations[orientation].append(runIndex) # Adding index to orientation list
        except FileNotFoundError as fileError:
            print(f"The file {fileError.filename} could not be opened for run {runIndex}!")
            print(f"ERROR: {fileError.strerror}")
            print(f"Continuing to the rest of the script: {sys.argv[0]}")

    # Printing orientation information to file:
    if reg.NUMBER_OF_POLYMER_REGIONS == 2: # Only printing orientation if there are two regions per polymer
        ratio = len(orientations["Parallel"]) / (len(orientations["Antiparallel"]) + 0.000001) # avoiding divide by 0 error with the 0.000001
        orientationFilePath = f"{GetFolder(sysPaths.CREATE_INITIAL_STATES, numberOfMonomers, architecture)}Analysis/orientation.txt"
        with open(orientationFilePath, "w") as file:
            for orient, runs in orientations.items():
                file.write(f"{orient}:\n{runs}\n")
            file.write(f"\nParallel to Antiparallel ratio = {len(orientations['Parallel'])}/{len(orientations['Antiparallel'])} = {ratio}\n")

def GetRegionCount(regionIndex: int) -> int:
    """Returns the number of monomers present in a particular region specified in region_config.py.
    The regionIndex is indexed from 0."""
    return (reg.REGION_BOUNDS[regionIndex][1] - reg.REGION_BOUNDS[regionIndex][0]) % numberOfMonomers + 1

def PlotPolymerCoMsSeparately(savePath: str, runIndex: int, total_CoM_List: list[ArrayLike], scaledTimeSteps: ArrayLike, axisLabel: str, showPlot: bool) -> None:
    """Plots the time series of the CoMs of the two polymers in the same plot and saves the figure.
    The total_CoM_List is a list of the numpy arrays containing the CoMs at each time step."""
    # Defining label to add to the y axis in case the distances are scaled with L:
    yLabelModifier = ""
    if segParam.RESCALE_LENGTHS:
        total_CoM_List = pt.RescaleArraysByLength(total_CoM_List, boxLength) # Rescaling lengths
        total_CoM_List, yAxisLabel = pt.ConvertMultipleArraysToScientificNotation(total_CoM_List)
        if len(yAxisLabel) != 0:
            yAxisLabel = f"({yAxisLabel})"
        yLabelModifier = rf"/ $L$ {yAxisLabel}"
    # Plotting: CoMs separately
    fig, ax = plt.subplots()
    for n in range(numberOfPolymers):
        ax.plot(scaledTimeSteps, total_CoM_List[n], label = f"Polymer {n+1}")
    if segParam.SHOW_TITLE:
        ax.set_title(f"Time series of CoMs of {numberOfPolymers} {architecture} polymers\n{numberOfMonomers} monomers each, Run {runIndex}")
    ax.set_xlabel(rf"Time ({axisLabel}$\tau_0$)")
    ax.set_ylabel(r"$z_{CoM}$ %s" % (yLabelModifier))
    leg = ax.legend()
    # Increasing legend line thickness:
    for line in leg.get_lines():
        line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)
    
    # saving figure:
    fig.savefig(savePath)
    plt.close(fig)
    print(f"Plotted Polymer CoM Time Series for {architecture} Run {runIndex}")


def PlotCoMDistance(savePath: str, runIndex: int, total_CoM_List: list[ArrayLike], scaledTimeSteps: ArrayLike, axisLabel: str, showPlot: bool) -> None:
    """Plots the CoM Distance after calculating it from the passed CoM list. It is plotted against the passed time steps on the passed axis.
    The time steps are expected to be rescaled by the LJ time.
    The plot is saved as the passed file path."""
    # Plotting: CoM Distance
    fig, ax = plt.subplots()
    # for n in range(numberOfPolymers):
    #     ax.plot(scaledTimeSteps, total_CoM_List[n], label = f"Polymer {n+1}")
    CoM_Distance = np.abs(total_CoM_List[0] - total_CoM_List[1])
    # Defining label to add to the y axis in case the distances are scaled with L:
    yLabelModifier = ""
    if segParam.RESCALE_LENGTHS:
        CoM_Distance = CoM_Distance / boxLength # Rescaling lengths
        CoM_Distance, yAxisLabel = pt.ConvertToScientificNotation(CoM_Distance)
        if len(yAxisLabel) != 0:
            yAxisLabel = f"({yAxisLabel})"
        yLabelModifier = rf" $L$ {yAxisLabel}"
    mean = np.mean(CoM_Distance)

    ax.plot(scaledTimeSteps, CoM_Distance, label = "Mean = %.3lf%s\nLast = %.3lf%s" % (mean, yLabelModifier, CoM_Distance[-1], yLabelModifier))
    if segParam.SHOW_TITLE:
        ax.set_title(f"Time series of CoM distance between {numberOfPolymers} {architecture} polymers\n{numberOfMonomers} monomers each, Run {runIndex}")
    else: # Printing a shorter title
        ax.set_title(segParam.GetPaperFiguresTitle(numberOfMonomers, architecture))
    ax.set_xlabel(rf"Time ({axisLabel}$\tau_0$)")
    if segParam.RESCALE_LENGTHS:
        yLabelModifier = rf"/ $L$ {yAxisLabel}"
    ax.set_ylabel(r"$\Delta z_{CoM}$ %s" % (yLabelModifier))
    leg = ax.legend()
    # Increasing legend line thickness:
    for line in leg.get_lines():
        line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)
    
    # saving figure:
    fig.subplots_adjust(bottom = 0.15)
    fig.savefig(savePath)
    plt.close(fig)
    print(f"Plotted Polymer CoM Distance for {architecture} Run {runIndex}")

def PlotPolymerCoMFromRegions(runIndex: int, showPlot: bool = False) -> None:
    """Plots the CoM time series of the two polymers after reading the region wise CoMs for a particular run"""
    folder = GetFolder(sysPaths.CREATE_INITIAL_STATES, numberOfMonomers, architecture)
    timeSteps, z_CoM_List = plotCoM.ReadRegionData(folder, runIndex)
    
    # Scaling timesteps by tau:
    scaledTimeSteps = timeSteps / TAU_0
    scaledTimeSteps, axisLabel = pt.ConvertToScientificNotation(scaledTimeSteps)

    # calculating the total polymer CoMs based on the number of monomers in each region
    regionsInSinglePolymer = reg.NUMBER_OF_REGIONS // numberOfPolymers # The number of regions in a single polymer
    total_CoM_List = [] # List of CoM arrays; one for each polymer
    for n in range(numberOfPolymers):
        polymerCoMs = np.zeros(len(z_CoM_List[0]))
        monomerCount = 0 # A variable to keep track of number of monomers in the regions
        for i in range(regionsInSinglePolymer):
            regionID = n * regionsInSinglePolymer + i
            polymerCoMs += GetRegionCount(i) * z_CoM_List[regionID]
            monomerCount += GetRegionCount(i)
        # Checking if the monomer counts add up to total:
        if numberOfMonomers != monomerCount:
            print(f"The monomers in regions ({monomerCount}) do not add up to {numberOfMonomers}!")
            sys.exit(1)
        total_CoM_List.append(polymerCoMs / numberOfMonomers) # normalizing the weighted sum

    # Plotting Time series of each polymer CoM separately in the same figure:
    # Plotting the CoM distance between two polymers:
    # Making the folder:
    CoMFolder = f"{folder}Analysis/Polymer_CoM/"
    Path(CoMFolder).mkdir(parents = True, exist_ok = True)
    PlotCoMDistance(f"{CoMFolder}CoM_Distance_r{runIndex}{segParam.FIG_EXT}", runIndex, total_CoM_List, scaledTimeSteps, axisLabel, showPlot)

def PlotPolymerRegionCoMForAllRuns(showPlots: bool = False):
    """Plots the polymer CoMs after reading the regional CoMs in the system for all runs"""
    for runIndex in range(1, numberOfRuns + 1):
        try:
            # timeSteps, z_com_list = ReadData(runIndex, readFileNamePrefix)
            # PlotData(timeSteps, z_com_list, showPlots)
            PlotPolymerCoMFromRegions(runIndex, showPlots)
        except FileNotFoundError as fileError:
            print(f"The file {fileError.filename} could not be opened for run {runIndex}!")
            print(f"ERROR: {fileError.strerror}")
            print(f"Continuing to the rest of the script: {sys.argv[0]}")

# Performing operations:

if __name__ == "__main__":
    SetArchitectureAndMonomers()
    SetConstants()
    # print(f"Length of Box: {boxLength}")
    runIndex = GetRunIndex()
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/big_bold.mplstyle")
    if runIndex == -1: # no run argument passed
        if reg.USE_REGIONS:
            PlotRegionTimeSeriesForAllRuns()
            PlotPolymerRegionCoMForAllRuns()
        else:
            PlotPolymerCoMForAllRuns()
    else:
        # try:
        #     readFileNamePrefix = "post_com"
        #     timeSteps, z_com_list = ReadData(runIndex, readFileNamePrefix)
        #     # PlotData(timeSteps, z_com_list, True)
        #     # Plotting post mixing CoM Distance:
        #     folder = sysPaths.CREATE_INITIAL_STATES
        #     PlotCoMDistance(timeSteps, z_com_list, folder)
        # except FileNotFoundError as fileError:
        #     print(f"One of {readFileNamePrefix} files could not be found!")
        #     print(f"ERROR: {fileError.filename} not found!")
        #     print("Continuing to the rest of the script")
        if reg.USE_REGIONS:
            ReadAndPlotRegionCoM(runIndex, True)
            PlotPolymerCoMFromRegions(runIndex, True)
        else:
            ReadAndPlotPolymerCoM(runIndex, True)
    print(f"Architecture: {architecture}")

# plotting CoM Distance:
# timeSteps, z_com_list = ReadData(runIndex)
# segTime = FindImprovedsegregationTIme(z_com_list)
# PlotCoMDistance(timeSteps, z_com_list, segTime)