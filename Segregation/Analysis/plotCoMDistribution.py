from io import TextIOWrapper
from os import times
import matplotlib as mpl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import pathlib
import Segregation_Parameters as segParam
from pathlib import Path

from typing import Tuple, List
from numpy.typing import ArrayLike
# importing the scripts that contains all the local paths
sys.path.append(f"../../Global_Scripts/System_File_Paths/") # Adding the path to the system file paths module to the system path
import system_file_paths as sysPaths
# Importing plotting tools:
sys.path.append(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/")
import PlottingTools as pt

# Importing the region Config file:
sys.path.append(f"{sysPaths.GLOBAL_SCRIPTS}Config_Files/")
import regions_config as reg
import plotMonomerDensity

# polymer information:
numberOfPolymers = 2
indexOfPolymer = 1
architecture = "Arc2"
runIndex = 1
boxLength = 25
initializationProcedure = "" # The name of the initialization procedure used to generate the initial state; this is used to read from the correct folder

numberOfRuns = 50
numberOfMonomers = 200 # each 
numberOfSteps = int(10 ** 7)
dataInterval = 100
numberOfBins = 65
acceptedModes = ["timeseries", "squared_timeseries", "runs_squared_timeseries", "squared_stddev", "distribution", "test_outliers"]
chosenMode = "timeseries" # can be "timeseries", "squared_timeseries", or "distribution"
numberOfMandatoryArguments = 3

def GetFolder(numberOfMonomers: int, architecture: str) -> str:
    """Returns the parent folder where all the segregation files lie"""
    prefix = sysPaths.GetFolder(sysPaths.SEGREGATION, numberOfMonomers, initializationProcedure)
    return f"{prefix}{architecture}/"

def GetPreviousAttemptsFolder(numberOfMonomers: int, architecture: str, simulationName: str) -> str:
    """Returns the path of the segregation times file in a particular backed-up smulation folder"""
    return f"{sysPaths.SEGREGATION}b{numberOfMonomers}/{simulationName}/{architecture}/"


def SetConstants():
    """Sets the global variables according to the case pertaining to the number of monomers used"""
    global boxLength
    global numberOfSteps
    global dataInterval

    diameter = pt.ReadDiameter(numberOfMonomers, architecture)
    boxLength = pt.GetAxisLength(numberOfMonomers, architecture, initializationProcedure, diameter)
    print(f"Diameter: {diameter}, Box length: {boxLength}")
    if numberOfMonomers == 200:
        numberOfSteps = int(4 * 10 ** 7)
        dataInterval = 1000
    elif numberOfMonomers == 500:
        numberOfSteps = int(10 ** 8)
        dataInterval = 10000

def SetArchitectureAndMonomers():
    """Reads the argument(s) passed while invoking the script from the command line and sets it to the architecture and numberOfMonomers"""
    global architecture
    global numberOfMonomers
    global initializationProcedure
    
    if len(sys.argv) < numberOfMandatoryArguments + 1:
        print("Not enough arguments specified while invoking the script! Please pass the number of monomers (integer), the architecture (string), and the initialization procedure (string) in the following format:")
        print("python /<path>/plotCoMDistribution.py <noOfMonomers> <architecture> <initializationProcedure>")
        print("An optional argument for the operation mode and/or the run index can be passed to plot the CoM data for only that run")
        print(f"Accepted operation modes: {acceptedModes}")
        quit() # terminating the script
    
    try: # checking if the number of monomers is valid:
        numberOfMonomers = int(sys.argv[1])
    except ValueError:
        print(f"The entered numberOfMonomers {sys.argv[1]} cannot be converted to a number! Please provide a valid number")
        quit()
    architecture = sys.argv[2]
    initializationProcedure = sys.argv[3]

def GetOptionalArguments():
    """
    Reads the last few optional arguments passed while invoking the script from the command line.
    Sets the chosen mode and returns the integer value for the run index.
    """
    global chosenMode
    if len(sys.argv) == numberOfMandatoryArguments + 1: # just script name, the architecture, and number of monomers
        chosenMode = "timeseries" # default mode
        return -1
    elif len(sys.argv) == numberOfMandatoryArguments + 2: # one optional argument passed
        arg = sys.argv[numberOfMandatoryArguments + 1]
        if arg.isdigit():
            return int(arg) # returns the run Index if the optional argument is an integer
        else: # the argument is a string
            if arg in acceptedModes:
                chosenMode = arg
                return -1
            else:
                print(f"The third argument was expected to be either an integer (run index) or one of the following strings: {acceptedModes}. But the one provided {arg} is neither.")
                sys.exit(1)
    else: # atleast two optional arguments passed
        # assuming the first optional argument is the mode and the second is the run index
        arg = sys.argv[numberOfMandatoryArguments + 1]
        if arg in acceptedModes:
            chosenMode = arg
        else:
            print(f"The first optional argument was expected to be one of the following strings: {acceptedModes}. But the one provided {arg} is not.")
            sys.exit(1)
        arg = sys.argv[numberOfMandatoryArguments + 2]
        if arg.isdigit():
            return int(arg) # returns the run Index if the optional argument is an integer
        else:
            print(f"The second optional argument was expected to be an integer for the run index, but the one provided {arg} cannot be cast to an integer")
            sys.exit(1)

def PlotTimeSeries(z_com: list, ax, indexOfPolymer: int, unitsAxisLabel: str = 'iterations', timeSteps = None):
    # plotting time series:
    if timeSteps is None:
        timeSteps = range(0, numberOfSteps, dataInterval)
    # Rescaling the CoM distances with box length
    yLabelModifier = ""
    if segParam.RESCALE_LENGTHS:
        z_com = z_com / boxLength
        yLabelModifier = r"/ $L$"

    ax.plot(timeSteps, z_com, label = f"Polymer {indexOfPolymer}")

    ax.set_xlabel(f"Time Steps ({unitsAxisLabel})")
    ax.set_ylabel(r"$z_{CoM}$ %s" % (yLabelModifier))
    if segParam.SHOW_TITLE:
        ax.set_title(f"Time Evolution of centre of mass\n {numberOfPolymers} {architecture} polymers, Run {runIndex}")
    else:
        ax.set_title(segParam.GetPaperFiguresTitle(numberOfMonomers, architecture))

def PlotCoMDistance(ax, timeSteps, z_com_list, unitsAxisLabel: str = "iterations", plotLabel: str = "") -> None:
    """Plots the CoM Distance on the passed axis.
    The unitsAxisLabel can be passed the units for the timesteps. Default is 'iterations'"""
    distance = np.abs(z_com_list[0] - z_com_list[1])
    # Rescaling the CoM distances with box length
    yLabelModifier = ""
    if segParam.RESCALE_LENGTHS:
        distance = distance / boxLength
        yLabelModifier = r"/ $L$"

    # Plotting:
    if segParam.USE_LOGLOG_PLOT:
        plottingFunction = ax.loglog
    elif segParam.USE_SEMILOG_PLOT:
        plottingFunction = ax.semilogx
    else:
        plottingFunction = ax.plot
    plottingFunction(timeSteps, distance, label = plotLabel)
    if segParam.SHOW_TITLE:
        ax.set_title(f"CoM Distance for two {architecture} polymers\nRun {runIndex}, {numberOfMonomers} monomers, {segParam.CRITERION_STRING}")
    # else:
    #     ax.set_title(segParam.GetPaperFiguresTitle(numberOfMonomers, architecture))
    ax.set_xlabel(rf"Time ({unitsAxisLabel})")
    ax.set_ylabel(r"$\Delta z_{COM}$ %s" % (yLabelModifier))


def PlotCoMDistanceAndSave(timeSteps, z_com_list, segTime: int | None = None, showPlot: bool = False):
    """Plots the distance between the centres of mass of the two polymers as a function of the timesteps and saves the figure"""
    fig, ax = plt.subplots()
    PlotCoMDistance(ax, timeSteps, z_com_list)
    if segTime != None:
        PlotVerticalLine(ax, segTime)

    if showPlot:
        plt.show(block = True)
    folder = GetFolder(numberOfMonomers, architecture)
    # fig.savefig(f"{folder}CoM_Distance_{architecture}_{numberOfMonomers}.png")
    fig.savefig(f"{folder}Analysis/CoM_Distance/CoM_Distance_r{runIndex}.png")

def PlotDistribution(z_com, ax, indexOfPolymer: int):
    # plotting:
    ax.hist(z_com, numberOfBins, alpha = 0.5, label = f"Polymer {indexOfPolymer}")
    ax.set_title(f"Histogram of Z coord of CoM\n {numberOfPolymers} {architecture} polymers, Run {runIndex}")
    ax.set_xlabel("z coordinate")
    ax.set_ylabel("Frequency")

def AddVerticalArrow(ax, xy: tuple[float, float], arrowStyle: str, color: str, text: str = "") -> None:
    """
    Adds a vertical arrow at the sepcified position xy with the specified properties.
    """
    arrowLength = 0.1 # Fraction of the plot height
    # Position of text:
    xytext = (xy[0], xy[1] + arrowLength * segParam.FIRST_THRESHOLD)
    ax.annotate(text, xy, xytext, arrowprops = dict(arrowstyle=arrowStyle, color = color))

def PlotVerticalLine(axes, time: int):
    axes.axvline(x = time, color = 'r', linestyle = '--', label = rf'Time of Segregation') # = {time} $\times 10^5\tau_0$')

def FindSegregationTime(z_com_list):
    # debugging:
    if type(z_com_list[0][0]) is type("string"):
        print("The read z CoM data was read as string! Some entries might not be float")
        

    # centre to centre distance:
    distances = np.absolute(z_com_list[0] - z_com_list[1])

    # criteria for segregation: the centre-centre distance must be greater than the threshold at atleast some time step
    threshold = boxLength/2
    segTime = -1
    for time, distance in zip(range(0, numberOfSteps, dataInterval), distances):
        if(distance > threshold):
            segTime = time
            break
    return segTime

def FindImprovedSegregationTime(z_com_list) -> tuple[int, int]:
    """The polymers are deemed to be segregated if the distance between their CoMs crosses a threshold f
    and the average distance remains above a threshold s for a subsequent interval of steps.
    Polymers are considered 'remixed' if their CoM distance falls below a threshold t, after which the process is restarted.
    This function returns a tuple of time steps:
    i) the step number of the first time at which f is crossed.
    ii) the step number at which the polymers are deemed to be segregated. """
    
    # TODO: Write the time at which the first threshold is reached for the first time to a file for each run
    # centre to centre distances:
    distances = np.absolute(z_com_list[0] - z_com_list[1])

    # flag to reflect segregated state:
    isSegregated = False

    firstThreshold = segParam.FIRST_THRESHOLD * boxLength # the threshold used for the first condition: checking if the CoM distance crosses this at any time
    secondThreshold = segParam.SECOND_THRESHOLD * boxLength # the threshold for the second condition: checking if the mean CoM distance is above this threshold; 0.5 * boxLength was too stringent
    thirdThreshold = segParam.THIRD_THRESHOLD * boxLength # the threshold to decide whether a segregated system gets mixed again; when the CoM distance falls below this threshold
    # assuming the CoM distance for a pair of segregated polymers would not go below this thresold just by fluctuations
    segTime = -1
    firstPassageTime = -1 # The time at which the CoM Distance crosses the first threshold for the first time
    # criterion for segregation: the initial point must be higher than threshold and the average of a subsequent interval too:
    intervalLength = int(numberOfSteps * segParam.INTERVAL_LENGTH) # in number of steps; half of total simulation run
    intervalLength = intervalLength // dataInterval # normalising to get the interval in terms of indices of the distance array
    for index, distance in zip(range(0, numberOfSteps // dataInterval), distances):

        if not isSegregated:
            if(distance > firstThreshold):
                # first condition met
                if firstPassageTime == -1: # Setting first passage time if not already set
                    firstPassageTime = index * dataInterval
                # the average of the subsequent interval:
                meanDistance = np.mean(distances[index: index + intervalLength])
                # debugging:
                # print(f"At time step {index * dataInterval}, mean = {meanDistance / boxlength} L")

                if meanDistance > secondThreshold:
                    # seconf condition met
                    segTime = index * dataInterval
                    isSegregated = True
                    # break # not breaking after segregation occurs; checking if it mixes back again
                else:
                    print(f"Segregation rejected at step {index * dataInterval} since average was below threshold")
                    continue
        else: # isSegregated
            #checking if it has mixed again
            if distance < thirdThreshold: # mixed again
                print(f"Segregated polymers mixed again! At step {index * dataInterval}")
                isSegregated = False
                segTime = -1 # resetting segTime
            else:
                continue
        
    return firstPassageTime, segTime

def ReadData(runIndex: int):
    """Reads the z coordinate of the CoM Data for both the polymers and returns the (timeSteps, double dimesnional list containing the read data) for a particular run"""
    global numberOfSteps
    folder = GetFolder(numberOfMonomers, architecture)
    z_com_list = []
    timeSteps = []
    isTimeStepsSet = False # a flag to indicate whether the time steps has been read and set already
    limitingSteps = numberOfSteps // dataInterval # the minium number of steps in the two CoM files; using this to make the shapes of both z_com arrays the same

    for indexOfPolymer in range(1, numberOfPolymers + 1):
        filePath = f"{folder}run{runIndex}/com{indexOfPolymer}.dat"
        df = pd.read_csv(filePath, sep = ' ')
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

    # Updating value of numberOfSteps to match the limitingSteps
    numberOfSteps = limitingSteps * dataInterval

    return timeSteps, z_com_list

def ScaleSegregationTimes(timeSteps: list[int], segTime: int) -> tuple[list[int], int, str]:
    """Rescales the passed timesteps and segregation time by Tau_0 and converts to scientific notation (if required).
    Returns a tuple containing:
    -Rescaled timesteps
    -Rescaled segregation time
    -Axis label or order of 20 factor in scientific notation"""
    # Coverting timesteps to units of TAU_0
    timeSteps = np.array(timeSteps) / segParam.TAU_0
    segTime = segTime / segParam.TAU_0
    scaledTimes = [timeSteps, np.array([segTime])] # Just putting the timesteps and the segTime in one array so they can be scaled topgether
    if not segParam.USE_LOGLOG_PLOT and not segParam.USE_SEMILOG_PLOT:
        scaledTimes, axisLabel = pt.ConvertMultipleArraysToScientificNotation(scaledTimes)
    else:
        axisLabel = ""
    timeSteps = scaledTimes[0]
    segTime = scaledTimes[1][0] # first element of the second array
    return (timeSteps, segTime, axisLabel)

def PlotData(timeSteps, z_com_list, firstPassageTime, segTime, showPlot: bool = False):
    """Takes the timesteps and the z_com data of the two polymers and plots one against the other. Plots a vertical red line at the time = segTime"""
    folder = GetFolder(numberOfMonomers, architecture)
    # fig1, ax1 = plt.subplots()
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/big_bold.mplstyle")
    fig2, ax2 = plt.subplots()

    didPolymersSegregate = True
    if segTime == -1:
        didPolymersSegregate = False
    # Coverting timesteps to units of TAU_0
    timeSteps, segTime, axisLabel = ScaleSegregationTimes(timeSteps, segTime)

    # average:
    # z_com = (z_com1 + z_com2)/2
    # PlotDistribution(z_com, ax1, indexOfPolymer)
    if segParam.SPLIT_COM_DISTANCE:
        for index in range(numberOfPolymers):
            PlotTimeSeries(z_com_list[index], ax2, index+1, rf"{axisLabel}$\tau_0$", timeSteps)
    else:
        PlotCoMDistance(ax2, timeSteps, z_com_list, rf"{axisLabel}$\tau_0$")

    # PlotCoMDistance(z_com_list[0], z_com_list[1], ax2, timeSteps)
    if didPolymersSegregate:
        PlotVerticalLine(ax2, segTime)
        # ax1.legend()
        # Plotting arrow for first passage time:
        if firstPassageTime != segTime:
            AddVerticalArrow(ax2, (firstPassageTime, segParam.FIRST_THRESHOLD), arrowStyle='->', color = 'gray', text = r'T_{F1}')
        ax2.legend(loc = (0.3, 0.7))

    if segParam.PLOT_INSET: # Plotting segregation trajectory comparisons of AOI_COMPARE as an inset figure
        # Creating the inset axis:
        mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/inset.mplstyle")
        axin = ax2.inset_axes(segParam.INSET_POSITION)
        PlotComparisonInSinglePlot(axin, runIndex)
        # Setting axis labels and legend:
        axin.set_ylabel("")
        axin.set_xlabel("")
        axin.set_xlim(left = -0.01, right = 0.4)
        axin.legend()

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig2, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)
    # fig1.savefig(f"{folder}Analysis/CoM_Distribution/CoMDistribution_r{runIndex}.png")
    # fig2.savefig(f"{folder}Analysis/TimeSeries/CoMTimeSeries_r{runIndex}.png")
    # fig2.subplots_adjust(bottom=0.15)
    if segParam.USE_LOGLOG_PLOT:
        plotLabel = "loglog_"
    elif segParam.USE_SEMILOG_PLOT:
        plotLabel = "semilog_"
    else:
        plotLabel = ""
    if segParam.SPLIT_COM_DISTANCE:
        fig2.savefig(f"{folder}Analysis/TimeSeries/{plotLabel}Polymer_CoM_r{runIndex}{segParam.FIG_EXT}")
    else:
        fig2.savefig(f"{folder}Analysis/TimeSeries/{plotLabel}CoMDistance_r{runIndex}{segParam.FIG_EXT}")
    
    #closing figures to save memory:
    # plt.close(fig1)
    plt.close(fig2)
    # print(f"Plotted CoM Time Series for {architecture} Run {runIndex}")
    print(f"Plotted CoM Distance for {architecture} Run {runIndex}")

def FindSegregationDirection(z_com_list) -> int:
    """Calculates the direction in which the segregation and occured for a particular run and returns an integer.
    Returns 0 if the 1st polymer ends up in the positive z-half of the cylinder.
    Returns 1 if the 1st polymer ends up in the negative z-half of the cylinder."""

    # Assuming that the polymers segregate during the run
    # Just comparing the z_com at the last timestep
    z_com_displacement = z_com_list[0][-1] - z_com_list[1][-1] # first polymer - second polymer
    if z_com_displacement > 0:
        return 0
    elif z_com_displacement < 0:
        return 1
    else:
        print("The com displacement is 0 at the last time step!")
        return 2 # should trigger out of bounds error

def WriteSegTimes(folder: str, segTimes: list[int], segDirections: list[int]) -> None:
    """Writes the passed segregation times to a file.
    Also writes the ratio of the number of runs segregating in different directions to another file.
    The argument folder is the directory path to the architecture folder.
    The argument segTimes is the list of segregation times for each run.
    The argument segDirections is the list of segregation directions for each run."""
    # array to count the segregation runs in one of the two directions:
    directionCounter = [0, 0] # 0th index is when 1st polymer is in +ve z half
    # and 1st index is when 1st polymer is in the negative z half

    # file to write the segregation times:
    segTimeFilePath = f"{folder}Analysis/segregationTimes_{segParam.CRITERION_STRING}.csv"
    segTimeFile = open(segTimeFilePath, "w")
    segTimeFile.write(f"Run Index, Segregation Time steps ({segParam.CRITERION_STRING})\n") # Header

    for runIndex in range(numberOfRuns):
        segTime = segTimes[runIndex]
        if(segTime == -1):
            print(f"Run {runIndex + 1} did not show segregation\n")
        else:
            segTimeFile.write(f"{runIndex + 1}, {segTime}\n")
            # print(f"{runIndex + 1}, {segTime}\n")
            try:
                directionCounter[segDirections[runIndex]] += 1
            except IndexError:
                print(f"For Run {runIndex}")
    segTimeFile.close()
    
    # Writing the ratio of segregated directions to a file:
    directionRatioFile = f"{folder}Analysis/direction_ratio.txt"
    with open(directionRatioFile, "w") as file:
        file.write("The ratio of 'right' segregations to total number of segregations:\n")
        file.write(f"{directionCounter[0]/sum(directionCounter)}\n")
        file.write("'Right' segregation: 1st polymer segregates to the positive z half\n")
        file.write("'Left' segregation: 1st polymer segregates to the negative z half")

def WriteFirstPassageTimes(folder: str, firstTimes: list[int]) -> None:
    """Writes the passed list of first passage times for each run to a new file under the folder directory
    The folder is directory path to the architecture folder."""
    # File to write first passage times:
    firstTimeFilePath = f"{folder}Analysis/firstPassageTimes_{segParam.CRITERION_STRING}.csv"
    with open(firstTimeFilePath, "w") as firstTimeFile:
        firstTimeFile.write(f"Run Index, First Passage Time steps ({segParam.CRITERION_STRING})\n") # Header
        for runIndex in range(numberOfRuns):
            firstPassageTime = firstTimes[runIndex]
            if firstPassageTime == -1:
                print(f"Run {runIndex + 1} never crossed the first threshold.")
                continue
            else:
                firstTimeFile.write(f"{runIndex + 1}, {firstPassageTime}\n")

def SegregateForAllRuns():
    global runIndex

    segTimes = [] # List to store segregation times
    segDirections = [] # List to store segregation directions: in which direction of the two each run segregated
    firstTimes = [] # List to store first passage times

    for runIndex in range(1, numberOfRuns + 1):
        timeSteps, z_com_list = ReadData(runIndex)
        firstPassageTime, segTime = FindImprovedSegregationTime(z_com_list)
        firstTimes.append(firstPassageTime)
        segTimes.append(segTime)
        segDirections.append(FindSegregationDirection(z_com_list))
        PlotData(timeSteps, z_com_list, firstPassageTime, segTime)

    # Writing times to files:
    folder = GetFolder(numberOfMonomers, architecture)
    WriteSegTimes(folder, segTimes, segDirections)
    WriteFirstPassageTimes(folder, firstTimes)

    print(f"Succesfully plotted the CoM Time Series for the {numberOfRuns} runs of {architecture}")

def PlotSquaredDeltaCOMDisplacementOnAxis(ax: mpl.axes.Axes, runIndex: int = -1, label: str = "", linestyle: str = "-", marker: str = '.', plotDiffusionGuide: bool = False, plotAllRuns: bool = False, skipOutlierRuns: bool = False, normalizeTime: bool = False, normalizeByMedianTime: bool = False) -> None:
    """
    Plots the squared delta COM displacement as a function of time i.e. (∆z_COM(t) - ∆z_COM(0))^2
    Args:
        - ax: The matplotlib axis on which the data is to be plotted
        - runIndex: The index of the run for which the squared displacement is to be plotted. If no runIndex is passed, the mean over all runs is plotted.
        - label: The label to be used for the plotted data
        - plotDiffusionGuide: Whether to plot a guide line for standard diffusion (z^2 = t) or not
        - plotAllRuns: Whether to plot all the individual runs in the background as light, grey lines or not
        - skipOutlierRuns: Whether to skip runs that contribute seg time outliers or not
        - normalizeTime: Whether to normalize the time by the segregation time or not for each individual run
        - normalizeByMedianTime: Whether to normalize the time by the median segregation time of all runs or not
    """
    if normalizeTime and normalizeByMedianTime:
        print("ERROR: Both normalizeTime and normalizeByMedianTime cannot be true at the same time!")
        sys.exit(1)
    timeStepsArray = []
    runDistances = []

    runIndices, segTimes = segParam.ReadSegregationTimes(numberOfMonomers, architecture, segParam.CRITERION_STRING)
    outlierIndices = segParam.FindOutlierIndices(segTimes)
    outlierRuns = [runIndices[i] for i in outlierIndices]
    if runIndex == -1: # averaging over all runs
        for i in range(1, numberOfRuns + 1):
            timeSteps, z_com_list = ReadData(i)
            timeStepsArray.append(timeSteps) # appending the timeSteps for each run to a list
            if segParam.COM_MODE == "Distance":
                distances = (z_com_list[0] - z_com_list[1]) ** 2 # (Delta z_COM(t))^2
            else:
                distances = (z_com_list[0] - z_com_list[1] - (z_com_list[0][0] - z_com_list[1][0])) ** 2 # Delta z_COM(t) - Delta z_COM(0)
            # Plotting each run individually to understand how the average comes about:
            runDistances.append(distances) # appending the distances for each run to a list
            if i == 1:
                allDistances = np.zeros(len(distances))
            allDistances += distances
        distance = allDistances / numberOfRuns
    else:
        if skipOutlierRuns and runIndex in outlierRuns:
            print(f"Skipping outlier {architecture} Run {runIndex} for Sq COM trajectory plot")
            return
        timeSteps, z_com_list = ReadData(runIndex)
        if segParam.COM_MODE == "Distance":
            distance = (z_com_list[0] - z_com_list[1]) ** 2 # (Delta z_COM(t))^2
        else:
            distance = (z_com_list[0] - z_com_list[1] - (z_com_list[0][0] - z_com_list[1][0])) ** 2 # Delta z_COM(t) - Delta z_COM(0)
    # Rescaling the CoM distances with box length
    yLabelModifier = r"($\sigma ^2$)"
    if segParam.RESCALE_LENGTHS:
        distance = distance / boxLength**2
        if runIndex == -1 and plotAllRuns:
            for i in range(numberOfRuns):
                runDistances[i] = runDistances[i] / boxLength**2
        yLabelModifier = r"/ $L^2$"
    
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/big_bold.mplstyle")
    # Scaling times:
    if normalizeTime:
        if runIndex == -1: # normalizing for all runs:
            for i in range(numberOfRuns):
                if (i+1) not in runIndices: # Did not segregate
                    timeStepsArray[i] = []
                else:
                    timeStepsArray[i] = timeStepsArray[i] / segTimes[i]
    elif normalizeByMedianTime:
        medianSegTime = np.median(segTimes)
        timeSteps = timeSteps / medianSegTime
    else:
        timeSteps = timeSteps / segParam.TAU_0
    if not segParam.USE_SEMILOG_PLOT and not segParam.USE_LOGLOG_PLOT:
        timeSteps, axisLabel = pt.ConvertToScientificNotation(timeSteps)
    else:
        axisLabel = ""

    # Plotting:
    if segParam.USE_LOGLOG_PLOT:
        plottingFunction = ax.loglog
    elif segParam.USE_SEMILOG_PLOT:
        plottingFunction = ax.semilogx
    else:
        plottingFunction = ax.plot
    # Plotting each run:
    if runIndex == -1 and plotAllRuns:
        for i in range(1, numberOfRuns + 1):
            if skipOutlierRuns and (i in outlierRuns):
                print(f"Skipping outlier {architecture} Run {i} for Sq COM trajectory plot")
                continue
            if normalizeTime and len(timeStepsArray[i-1]) == 0:
                print(f"Run {i} did not show segregation; skipping for normalization of times")
                continue
            if normalizeTime:
                plottingFunction(timeStepsArray[i-1], runDistances[i-1], color = 'gray', alpha = 0.5, rasterized = True)
            else:
                plottingFunction(timeSteps, runDistances[i-1], color = 'gray', alpha = 0.5, rasterized = True) # Plotting unnormalized times
    if not normalizeTime:
        plottingFunction(timeSteps, distance, marker = marker, label = label, linestyle = linestyle, rasterized = True) # Does not make sense to plot mean for normalized runs
 
    if normalizeTime:
        ax.set_xlabel(r"Time / $T_{seg}$%s" % (axisLabel))
    elif normalizeByMedianTime:
        ax.set_xlabel(r"Time / Median $T_{seg}$%s" % (axisLabel))
    else:
        ax.set_xlabel(rf"Time ({axisLabel}$\tau_0$)")
    ax.set_ylabel(r"$(\Delta z_{COM})^2$ %s" % (yLabelModifier))
    # Plotting a guide line for standard diffusion for infinite cylinder situations
    if plotDiffusionGuide:
        prefactor, exponent = (0.1, 1)
        diffusionGuide = prefactor * timeSteps ** exponent
        ax.plot(timeSteps, diffusionGuide, linestyle = '--', color = 'red', label = f"Diffusion Guide: Exponent = {exponent}")
        ax.legend()

def PlotSqDeltaCOMDisplacement(runIndex: int = -1, showPlot: bool = False) -> None:
    """
    Plots the squared delta CoM displacement.
    The runIndex can be passed to compute the squared distance for that run; otherwise a mean over all runs is plotted.
    """
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/big_bold.mplstyle")
    fig, ax = plt.subplots()
    plotDiffusionGuide = False
    if sysPaths.IsCylinderInfinite() and segParam.USE_LOGLOG_PLOT and not segParam.NORMALIZE_TIME:
        plotDiffusionGuide = True

    PlotSquaredDeltaCOMDisplacementOnAxis(ax, runIndex, plotDiffusionGuide=plotDiffusionGuide, plotAllRuns=segParam.PLOT_ALL_RUNS, skipOutlierRuns=segParam.SKIP_OUTLIER_RUNS, normalizeTime=segParam.NORMALIZE_TIME)

    # Run label:
    if runIndex == -1:
        outlierRuns = 0
        outliersLabel = ""
        if segParam.SKIP_OUTLIER_RUNS:
            runIndices, segTimes = segParam.ReadSegregationTimes(numberOfMonomers, architecture, segParam.CRITERION_STRING)
            outlierIndices = segParam.FindOutlierIndices(segTimes)
            outlierRuns = len(outlierIndices)
            outliersLabel = " (skipped outliers)"
        runIndexLabel = f"Mean over {numberOfRuns-outlierRuns} runs{outliersLabel}"
    else:
        runIndexLabel = f"Run {runIndex}"
    if segParam.SHOW_TITLE:
        ax.set_title(f"Squared ∆CoM {segParam.COM_MODE} for two {architecture} polymers\n{runIndexLabel}, {numberOfMonomers} monomers, {segParam.CRITERION_STRING}")
    else:
        ax.text(0.1, 0.9, "%s\n%s\n%s" % (segParam.GetAliasArchitecture(architecture), rf"$N = {numberOfMonomers}$", segParam.GetAliasSimulation(sysPaths.SPECIAL_SIMULATION)), transform = ax.transAxes, horizontalalignment='left', verticalalignment='center')
    # ax.set_xlim(left = -0.01, right= 0.80)
    if showPlot:
        plt.show(block = True)

    # Saving figure:
    folder = GetFolder(numberOfMonomers, architecture)
    if runIndex == -1:
        saveRunLabel = 'all'
    else:
        saveRunLabel = f"r{runIndex}"
    fig.savefig(f"{folder}Analysis/TimeSeries/CoM_Sq_{segParam.COM_MODE}_{saveRunLabel}{segParam.FIG_EXT}")
    # Success message:
    print(f"Successfully plotted squared delta CoM {segParam.COM_MODE} for {runIndexLabel} of {architecture}")

def PlotSqDeltaCOMArcComparison(arcList: list[str], showPlot: bool = False, runIndex: int = -1, saveArcLabel: str = "arc_comp") -> None:
    """
    Plots the squared delta CoM Displacement for different architectures in a subplots figure.
    The architectures are specified in the arcList argument.
    Args:
        - arcList: The list of architectures to be compared.
        - showPlot: Whether to show the plot or not.
        - runIndex: The index of the run for which the squared delta CoM displacement is to be plotted. If none is passed, it plots an average over all runs
        - saveRunLabel: The label to be used while saving the figure.
    """
    global architecture
    folder = GetFolder(numberOfMonomers, architecture)
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/big_bold.mplstyle")
    fig, ax = plt.subplots()
    for i in range(len(arcList)):
        architecture = arcList[i]
        if i == 0 and sysPaths.IsCylinderInfinite() and segParam.USE_LOGLOG_PLOT and not segParam.NORMALIZE_TIME:
            plotDiffusionGuide = True
        else:
            plotDiffusionGuide = False
        PlotSquaredDeltaCOMDisplacementOnAxis(ax, runIndex, label = arcList[i], plotDiffusionGuide=plotDiffusionGuide, plotAllRuns=False, skipOutlierRuns=segParam.SKIP_OUTLIER_RUNS, normalizeTime=False, normalizeByMedianTime=segParam.NORMALIZE_TIME)
    ax.legend()
    # Run label:
    if runIndex == -1:
        runIndexLabel = f"Mean over {numberOfRuns} runs"
    else:
        runIndexLabel = f"Run {runIndex}"
    if segParam.SHOW_TITLE:
        ax.set_title(f"Squared ∆CoM {segParam.COM_MODE} Comparison\n{numberOfMonomers} monomers, {runIndexLabel}")
    else:
        ax.text(0.1, 0.9, "%s\n%s" % (rf"$N = {numberOfMonomers}$", segParam.GetAliasSimulation(sysPaths.SPECIAL_SIMULATION)), transform = ax.transAxes, horizontalalignment='left', verticalalignment='center')

    if showPlot:
        plt.show()
    if runIndex == -1:
        saveRunLabel = 'all'
    else:
        saveRunLabel = f"r{runIndex}"
    if segParam.NORMALIZE_TIME:
        # Label for normalized time:
        saveArcLabel = f"{saveArcLabel}_normTime"
        # Setting axis limit
        # if sysPaths.IsCylinderInfinite() and segParam.USE_LOGLOG_PLOT:
        #     ax.set_ylim(bottom = 0.1)
            
    fig.savefig(f"{folder}Analysis/TimeSeries/CoM_Sq_{segParam.COM_MODE}_{saveArcLabel}_{saveRunLabel}{segParam.FIG_EXT}")
    plt.close(fig)
    print(f"Successfully plotted squared delta CoM {segParam.COM_MODE} for multiple architectures: {arcList}")

def PlotSqDeltaCOMSizeComparison(sizeList: list[int], arcList: list[str], runIndex: int = -1, showPlot: bool = False, saveSizeLabel: str = "size_comp") -> None:
    """
    Plots the squared delta CoM Displacement for different polymer sizes (numbers of monomers) in a subplots figure.
    The sizes are specified in the sizeList argument.
    The architectures are specified in the arcList argument.
    Args:
        - sizeList: The list of polymer sizes (number of monomers) to be compared.
        - arcList: The list of architectures to be compared.
        - showPlot: Whether to show the plot or not.
        - runIndex: The index of the run for which the squared delta CoM displacement is to be plotted. If none is passed, it plots an average over all runs
        - saveRunLabel: The label to be used while saving the figure.
    """
    global architecture, numberOfMonomers, diameter, boxLength
    segCriteria = ["f045_s040_t00", "f048_s043_t00"] # Segregation criteria for different polymer sizes
    folder = GetFolder(numberOfMonomers, architecture)
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/big_bold.mplstyle")
    fig, ax = plt.subplots()
    # Darker gradients for better visibility on white backgrounds
    blue_gradient = ["#003366", "#004080", "#0059b3", "#0073e6", "#3399ff"]
    red_gradient = ["#660000", "#800000", "#b30000", "#e60000", "#ff3333"]
    colourWheels = [blue_gradient, red_gradient]
    for i in range(len(sizeList)):
        numberOfMonomers = sizeList[i]
        diameter = pt.ReadDiameter(numberOfMonomers, architecture)
        boxLength = pt.GetAxisLength(numberOfMonomers, architecture, initializationProcedure, diameter)
        segParam.CRITERION_STRING = segCriteria[i % len(segCriteria)] # cycling through the criteria
        ax.set_prop_cycle(color = colourWheels[i % len(colourWheels)])
        for j in range(len(arcList)):
            architecture = arcList[j]
            label = ""
            if j == 0:
                label = f"{sizeList[i]}" # Only labelling first of the architectures
            PlotSquaredDeltaCOMDisplacementOnAxis(ax, runIndex, label = label, plotDiffusionGuide=False, plotAllRuns=False, skipOutlierRuns=segParam.SKIP_OUTLIER_RUNS, normalizeTime=False, normalizeByMedianTime=segParam.NORMALIZE_TIME)
    ax.legend()
    # Run label:
    if runIndex == -1:
        runIndexLabel = f"Mean over {numberOfRuns} runs"
    else:
        runIndexLabel = f"Run {runIndex}"
    if segParam.SHOW_TITLE:
        ax.set_title(f"Squared ∆CoM {segParam.COM_MODE} Comparison\n{runIndexLabel}")
    else:
        ax.text(0.1, 0.9, segParam.GetAliasSimulation(sysPaths.SPECIAL_SIMULATION), transform = ax.transAxes, horizontalalignment='left', verticalalignment='center')

    if showPlot:
        plt.show(block = True)
    if runIndex == -1:
        saveRunLabel = 'all'
    else:
        saveRunLabel = f"r{runIndex}"
    if segParam.NORMALIZE_TIME:
        # Label for normalized time:
        saveSizeLabel = f"{saveSizeLabel}_normTime"

    # Saving figure:
    fig.savefig(f"{folder}Analysis/TimeSeries/CoM_Sq_{segParam.COM_MODE}_{saveSizeLabel}_{saveRunLabel}{segParam.FIG_EXT}")
    plt.close(fig)
    print(f"Successfully plotted squared delta CoM {segParam.COM_MODE} for multiple architectures: {arcList} and multiple size: {sizeList}")


def PlotSquaredCOMRunComparison(runIndices: list[int], showPlot: bool = False, groupName: str = "", saveRunLabel: str = "run_comp") -> None:
    """
    Plots the squared delta CoM Displacement for all runs individually on the passed axis.
    Args:
        - runIndices: The list of run indices for which the squared delta CoM displacement is to be plotted.
        - showPlot: Whether to show the plot or not.
        - groupName: The name of the group of runs being plotted. This is just for labelling purposes.
        - saveRunLabel: The label to be used while saving the figure.
    """
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/big_bold.mplstyle")
    fig, ax = plt.subplots()
    for runIndex in runIndices:
        PlotSquaredDeltaCOMDisplacementOnAxis(ax, runIndex)
    if  segParam.SHOW_TITLE:
        ax.set_title(f"Squared ∆CoM {segParam.COM_MODE} for two {architecture} polymers\n{len(runIndices)}{groupName} runs, {numberOfMonomers} monomers, {segParam.CRITERION_STRING}")
    else:
        ax.text(0.05, 0.8, "%s\n$N=%i$\n%s" % (segParam.GetAliasArchitecture(architecture), numberOfMonomers, segParam.GetAliasSimulation(sysPaths.SPECIAL_SIMULATION)), verticalalignment='center', horizontalalignment='left', transform=ax.transAxes)

    if segParam.USE_MARGINS:
        fig.subplots_adjust(left=segParam.PLOT_MARGINS['left'], right=segParam.PLOT_MARGINS['right'], top=segParam.PLOT_MARGINS['top'], bottom=segParam.PLOT_MARGINS['bottom'])

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)

    # Saving figure:
    folder = GetFolder(numberOfMonomers, architecture)
    fig.savefig(f"{folder}Analysis/TimeSeries/CoM_Sq_{segParam.COM_MODE}_{saveRunLabel}{segParam.FIG_EXT}")
    # Success message:
    print(f"Successfully plotted squared delta CoM {segParam.COM_MODE} for {runIndices} of {architecture}")

def PlotSqDeltaCOMGroupsComparison(showPlot: bool = False) -> None:
    """
    Plots the squared delta CoM Displacement for a sample of runs in different groups in a subplots figure.
    The groups are specified in the segParam.GROUP_INDICES variable in Segregation_Parameters.py
    """
    folder = GetFolder(numberOfMonomers, architecture)
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/big_bold.mplstyle")
    fig, ax = plt.subplots()
    segParam.UpdateGroupIndices(set(range(1, numberOfRuns + 1)))
    runLimit = 5
    # Gradient of 5 blue shades and 5 red shades (hex values)
    blue_gradient = ["#0077FF", "#3399FF", "#66B2FF", "#99CCFF", "#CCE5FF"]
    # red_gradient = ["#FF3300", "#FF6633", "#FF9966", "#FFB299", "#FFD6CC"]
    # blue_gradient = ["#003366", "#004080", "#0059b3", "#0073e6", "#3399ff"]
    red_gradient = ["#660000", "#800000", "#b30000", "#e60000", "#ff3333"]
    colourWheels = [blue_gradient, red_gradient]
    markers = [None, '.'] 
    linestyles = [':', '-' ] # ['..', '-']
    # TODO: Try different symbols/linestyles
    for i in range(len(segParam.GROUP_INDICES)):
        # Setting the colour cycle for the group:
        ax.set_prop_cycle(color = colourWheels[i % len(colourWheels)])
        for j in range(min(len(segParam.GROUP_INDICES[i]), runLimit)):
            if j == 0:
                groupLabel = segParam.GROUP_LABELS[i]
            else:
                groupLabel = ""
            PlotSquaredDeltaCOMDisplacementOnAxis(ax, runIndex = segParam.GROUP_INDICES[i][j], label = groupLabel, linestyle = linestyles[i % len(linestyles)], marker = markers[i % len(markers)], plotDiffusionGuide=False, plotAllRuns=False, skipOutlierRuns=False, normalizeTime=False)
    leg = ax.legend(loc = 'lower right')
    # Increasing legend line thickness:
    for line in leg.get_lines():
        line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())
        line.set_markersize(pt.LEGEND_LINEWIDTH * line.get_markersize())
    if segParam.SHOW_TITLE:
        ax.set_title(f"Squared ∆CoM {segParam.COM_MODE} Comparison for two {architecture} polymers\n{numberOfMonomers} monomers,  {len(segParam.GROUP_INDICES)} groups of runs")
    else:
        ax.text(0.05, 0.8, "%s\n%s\n%s" % (segParam.GetAliasArchitecture(architecture), rf"$N = {numberOfMonomers}$", segParam.GetAliasSimulation(sysPaths.SPECIAL_SIMULATION)), transform = ax.transAxes, verticalalignment='center', horizontalalignment='left')

    if segParam.USE_MARGINS:
        fig.subplots_adjust(left=segParam.PLOT_MARGINS['left'], right=segParam.PLOT_MARGINS['right'], top=segParam.PLOT_MARGINS['top'], bottom=segParam.PLOT_MARGINS['bottom'])

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show()
    fig.savefig(f"{folder}Analysis/TimeSeries/CoM_Sq_{segParam.COM_MODE}_groups{segParam.FIG_EXT}")
    plt.close(fig)
    print(f"Successfully plotted squared delta CoM {segParam.COM_MODE} for multiple groups of runs: {segParam.GROUP_INDICES}")

def CalculateSqCOMStdDev(normalize: bool = False, skipOutlierRuns: bool = False) -> tuple[ArrayLike, ArrayLike]:
    """
    Calculates the standard deviation of the squared delta CoM displacement over all runs for each timestep of the trajectories.
    Returns the timeSteps and the standard deviation as a tuple.
    The timesteps are in units of the simulation.
    Args:
        - normalize: Whether to normalize the standard deviation by the mean squared distance or not. 
                    If false, the standard deviation is in units of length in the simulation.
                    If true, the standard deviation is unitless, and standardized.
        - skipOutlierRuns: Whether to skip outlier runs while calculating the standard deviation or not.
    """
    # Reading segregation times and finding outlier runs:
    runIndices, segTimes = segParam.ReadSegregationTimes(numberOfMonomers, architecture)
    outlierIndices = segParam.FindOutlierIndices(segTimes)
    outlierRuns = [runIndices[i] for i in outlierIndices]
    allDistances = []
    for i in range(1, numberOfRuns + 1):
        if skipOutlierRuns and i in outlierRuns:
            print(f"Skipping outlier {architecture} Run {i} for stdDev calculation")
            continue
        timeSteps, z_com_list = ReadData(i)
        if segParam.COM_MODE == "Distance":
            distances = (z_com_list[0] - z_com_list[1]) ** 2 # (Delta z_COM(t))^2
        elif segParam.COM_MODE == "Displacement":
            distances = (z_com_list[0] - z_com_list[1] - (z_com_list[0][0] - z_com_list[1][0])) ** 2 # (Delta z_COM(t))^2
        allDistances.append(distances)
    allDistances = np.array(allDistances)
    stdDev = np.std(allDistances, axis = 0) # Axis 0: calculates std across all runs for a particular timestep
    meanDistances = np.mean(allDistances, axis = 0)
    if normalize:
        stdDev = stdDev / meanDistances

    return timeSteps, stdDev

def PlotSqCOMStdDev(normalize: bool = False, showPlot: bool = False) -> None:
    """
    Plots the standard deviation of the squared delta CoM displacement over all runs.
    Args:
        - normalize: Whether to normalize the standard deviation by the mean squared distance or not. 
                    If false, the standard deviation is in units of length in the simulation.
                    If true, the standard deviation is unitless, and standardized.
        - showPlot: Whether to show the plot or not.
    """
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/big_bold.mplstyle")
    fig, ax = plt.subplots()
    timeSteps, stdDev = CalculateSqCOMStdDev(normalize, segParam.SKIP_OUTLIER_RUNS) # In units of the simulation

    # Rescaling the CoM distances with box length
    yLabelModifier = ""
    if normalize:
        yLabelModifier = r"/ $<(\Delta z_{COM})^2>$"
    elif segParam.RESCALE_LENGTHS:
        stdDev = stdDev / boxLength**2
        yLabelModifier = r"/ $L^2$"

    # Scaling times:
    timeSteps = timeSteps / segParam.TAU_0
    if not segParam.USE_SEMILOG_PLOT and not segParam.USE_LOGLOG_PLOT:
        timeSteps, axisLabel = pt.ConvertToScientificNotation(timeSteps)
    else:
        axisLabel = ""

    # Plotting:
    if segParam.USE_LOGLOG_PLOT:
        plottingFunction = ax.loglog
    elif segParam.USE_SEMILOG_PLOT:
        plottingFunction = ax.semilogx
    else:
        plottingFunction = ax.plot
    plottingFunction(timeSteps, stdDev, marker = '.')
 
    ax.set_xlabel(rf"Time ({axisLabel}$\tau_0$)")
    ax.set_ylabel(r"$\sigma \left((\Delta z_{COM})^2 \right)$ %s" % (yLabelModifier))
    if segParam.SKIP_OUTLIER_RUNS:
        noOfOutliers = len(segParam.FindOutlierIndices(segParam.ReadSegregationTimes(numberOfMonomers, architecture)[1]))
        outliersSkippedLabel = " (outliers skipped)"
    else:
        noOfOutliers = 0
        outliersSkippedLabel = ""
    if segParam.SHOW_TITLE:
        ax.set_title(f"Standard Deviation of Squared ∆CoM {segParam.COM_MODE} for two {architecture} polymers\nMean over {numberOfRuns-noOfOutliers} runs{outliersSkippedLabel}, {numberOfMonomers} monomers, {segParam.CRITERION_STRING}")
    else:
        ax.text(0.1, 0.9, "%s\n%s\n%s" % (segParam.GetAliasArchitecture(architecture), rf"$N = {numberOfMonomers}$", segParam.GetAliasSimulation(sysPaths.SPECIAL_SIMULATION)), horizontalalignment='left', verticalalignment='center')

    if segParam.USE_MARGINS:
        fig.subplots_adjust(left=segParam.PLOT_MARGINS['left'], right=segParam.PLOT_MARGINS['right'], top=segParam.PLOT_MARGINS['top'], bottom=segParam.PLOT_MARGINS['bottom'])

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)

    # Saving figure:
    folder = GetFolder(numberOfMonomers, architecture)
    modifierLabel = ""
    if normalize:
        modifierLabel += "norm_"
    if segParam.SKIP_OUTLIER_RUNS:
        modifierLabel += "skip_"
    
    fig.savefig(f"{folder}Analysis/TimeSeries/CoM_Sq_{segParam.COM_MODE}_{modifierLabel}stdDev{segParam.FIG_EXT}")
    # Success message:
    print(f"Successfully plotted {modifierLabel}standard deviation of squared delta CoM {segParam.COM_MODE} for mean over all runs of {architecture}")

def PlotSqCOMStdDevCompOnAxis(ax: mpl.axes.Axes, arcList: list[str], normalize: bool = False, label: str = "", skipOutlierRuns: bool = False, plotMedianTime: bool = False, marker: str|None = None, fillstyle: str = "full", plottingInterval: int = 1) -> None:
    """
    Plots the standard deviation of the squared delta CoM displacement for different architectures in a subplots figure.
    The architectures are specified in the arcList argument.
    Args:
        - ax: The matplotlib axis on which the data is to be plotted
        - arcList: The list of architectures to be compared.
        - normalize: Whether to normalize the standard deviation by the mean squared distance or not. 
                    If false, the standard deviation is in units of length in the simulation.
                    If true, the standard deviation is unitless, and standardized.
        - label: The label to be used for the plotted data
        - skipOutlierRuns: Whether to skip outlier runs while calculating the standard deviation or not.
        - plotMedianTime: Whether to plot a vertical line at the median segregation time for each architecture or not.
        - marker: The marker to be used for the plotted data. If None, different markers for each architecture is used
        - fillstyle: The fillstyle to be used for the markers. Default is 'full'.
        - plottingInterval: The interval at which to plot the data points. Default is 1 (i.e. plot all points).
    """
    global architecture
    # Plotting style:
    if segParam.USE_LOGLOG_PLOT:
        plottingFunction = ax.loglog
    elif segParam.USE_SEMILOG_PLOT:
        plottingFunction = ax.semilogx
    else:
        plottingFunction = ax.plot
    for i in range(len(arcList)):
        architecture = arcList[i]
        timeSteps, stdDev = CalculateSqCOMStdDev(normalize, segParam.SKIP_OUTLIER_RUNS) # In units of the simulation
        # Scaling times:
        timeSteps = timeSteps / segParam.TAU_0
        if not segParam.USE_SEMILOG_PLOT and not segParam.USE_LOGLOG_PLOT:
            timeSteps, axisLabel = pt.ConvertToScientificNotation(timeSteps)
        else:
            axisLabel = ""
        # Rescaling the CoM distances with box length
        if not normalize and segParam.RESCALE_LENGTHS:
            stdDev = stdDev / boxLength**2
        # Get the color from the current color cycle for this architecture
        color = ax._get_lines.get_next_color()
        if len(label) != 0:
            if i != 0:
                plotLabel = ""
            else:
                plotLabel = label
        else:
            plotLabel = segParam.GetAliasArchitecture(arcList[i])
        # Setting marker:
        if marker is not None:
            plotMarker = marker
        else:
            plotMarker = segParam.MARKER_WHEEL[i]
        # print("Debugging: Plotting %s with marker %s and color %s. Length of array = %i" % (plotLabel, plotMarker, color, len(timeSteps[::plottingInterval])))
        plottingFunction(timeSteps[::plottingInterval], stdDev[::plottingInterval], marker = plotMarker, label = plotLabel, color=color, fillstyle = fillstyle, rasterized = True)
        # Plotting median segregation time as an indication for plot peak, using the same color
        if plotMedianTime:
            runIndices, segTimes = segParam.ReadSegregationTimes(numberOfMonomers, architecture, segParam.CRITERION_STRING)
            medianSegTime = np.median(segTimes)
            medianSegTime = medianSegTime / segParam.TAU_0
            ax.axvline(x = medianSegTime, linestyle = '--', linewidth = 2, color=color)
        print(f"Plotted Std. Dev. of Squared COM Distance for {sysPaths.SPECIAL_SIMULATION} {arcList[i]}")
    ax.set_xlabel(rf"Time ({axisLabel}$\tau_0$)")
    # Success message:

def PlotSqCOMStdDevArcComparison(arcList: list[str], normalize: bool = False, showPlot: bool = False) -> None:
    """
    Plots the standard deviation of the squared delta CoM displacement for different architectures in a subplots figure.
    The architectures are specified in the arcList argument.
    Args:
        - arcList: The list of architectures to be compared.
        - normalize: Whether to normalize the standard deviation by the mean squared distance or not. 
                    If false, the standard deviation is in units of length in the simulation.
                    If true, the standard deviation is unitless, and standardized.
        - showPlot: Whether to show the plot or not.
    """
    global architecture
    folder = GetFolder(numberOfMonomers, architecture)
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/big_bold.mplstyle")
    mpl.rcParams['axes.formatter.limits'] = (-3, 3) # To avoid scientific notation on axes
    fig, ax = plt.subplots()
    PlotSqCOMStdDevCompOnAxis(ax, arcList, normalize, skipOutlierRuns=segParam.SKIP_OUTLIER_RUNS, plotMedianTime = True, marker = '.')
    
    leg = ax.legend()
    # Increasing legend line thickness:
    for line in leg.get_lines():
        line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())
        line.set_markersize(pt.LEGEND_LINEWIDTH * line.get_markersize())
    if segParam.SKIP_OUTLIER_RUNS:
        outliersSkippedLabel = " (outliers skipped)"
    else:
        outliersSkippedLabel = ""
    if segParam.SHOW_TITLE:
        ax.set_title(f"Standard Deviation of Squared ∆CoM {segParam.COM_MODE} Comparison\n{numberOfMonomers} monomers, over all runs{outliersSkippedLabel}")
    else:
        ax.text(0.98, 0.9, "%s\n%s" % (rf"$N = {numberOfMonomers}$", segParam.GetAliasSimulation(sysPaths.SPECIAL_SIMULATION)), transform = ax.transAxes, horizontalalignment='right', verticalalignment='center')
    yLabelModifier = ""
    if normalize:
        yLabelModifier = r"/ $<(\Delta z_{COM})^2>$"
    elif segParam.RESCALE_LENGTHS:
        yLabelModifier = r"/ $L^2$"
    ax.set_ylabel(r"$\sigma \left((\Delta z_{COM})^2 \right)$ %s" % (yLabelModifier))

    if segParam.USE_MARGINS:
        fig.subplots_adjust(left=segParam.PLOT_MARGINS['left'], right=segParam.PLOT_MARGINS['right'], top=segParam.PLOT_MARGINS['top'], bottom=segParam.PLOT_MARGINS['bottom'])

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show()
    modifierLabel = ""
    if normalize:
        modifierLabel += "norm_"
    if segParam.SKIP_OUTLIER_RUNS:
        modifierLabel += "skip_"
    fig.savefig(f"{folder}Analysis/TimeSeries/CoM_Sq_{segParam.COM_MODE}_{modifierLabel}stdDev_arcComp{segParam.FIG_EXT}")
    plt.close(fig)
    print(f"Successfully plotted {modifierLabel}standard deviation of squared delta CoM {segParam.COM_MODE} for multiple architectures: {arcList}")

def PlotSqCOMStdDevInitComparison(initList: list[str], arcList: list[str], normalize: bool, showPlot: bool) -> None:
    """
    Plots the Sqaured COM Standard Deviation for multiple initialization procedures and multiple architectures.
    """
    global architecture
    saveFolder = GetFolder(numberOfMonomers, architecture)
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/big_bold.mplstyle")
    # Define a colour wheel with 5 contrasting and visible colours for a white background
    fig, ax = plt.subplots()
    fillStyles = ["full", "none"]
    colourWheel = mpl.rcParams['axes.prop_cycle'].by_key()['color']
    for i in range(len(initList)):
        # colourWheel = ["blue", "orange", "green", "red", "black"]  # blue, orange, green, red, black
        # Reset the color cycle to start from the first color again
        ax.set_prop_cycle(color=colourWheel)
        # ax._get_lines.prop_cycler = iter(ax._get_lines.prop_cycler)
        sysPaths.SPECIAL_SIMULATION = initList[i]
        folder = GetFolder(numberOfMonomers, architecture)
        PlotSqCOMStdDevCompOnAxis(ax, arcList, normalize = normalize, label = segParam.GetAliasSimulation(initList[i]), skipOutlierRuns = segParam.SKIP_OUTLIER_RUNS, marker = segParam.MARKER_WHEEL[i], fillstyle=fillStyles[i % len(fillStyles)], plottingInterval = 1, plotMedianTime=True)
    # Making custom legend:
    leg = ax.legend()
    # Get existing handles and labels
    handles, labels = leg.legend_handles, [t.get_text() for t in leg.get_texts()]

    # Creating proxy legend handles:
    proxies = []
    for i in range(len(arcList)):
        proxy = mpl.lines.Line2D([], [], color=colourWheel[i], linestyle='-', label=segParam.GetAliasArchitecture(arcList[i]))
        proxies.append(proxy)

    # Append new ones
    handles += proxies
    labels  += [p.get_label() for p in proxies]

    # Update legend
    ax.legend(handles, labels)

    if segParam.SKIP_OUTLIER_RUNS:
        outliersSkippedLabel = " (outliers skipped)"
    else:
        outliersSkippedLabel = ""
    if segParam.SHOW_TITLE:
        ax.set_title(f"Standard Deviation of Squared ∆CoM {segParam.COM_MODE} Comparison\n{numberOfMonomers} monomers, over all runs{outliersSkippedLabel}")
    else:
        ax.text(0.1, 0.9, "%s" % (rf"$N = {numberOfMonomers}$"), transform=ax.transAxes, horizontalalignment='left', verticalalignment='center')
    yLabelModifier = ""
    if normalize:
        yLabelModifier = r"/ $<(\Delta z_{COM})^2>$"
    elif segParam.RESCALE_LENGTHS:
        yLabelModifier = r"/ $L^2$"
    ax.set_ylabel(r"$\sigma \left((\Delta z_{COM})^2 \right)$ %s" % (yLabelModifier))

    if segParam.USE_MARGINS:
        fig.subplots_adjust(left=segParam.PLOT_MARGINS['left'], right=segParam.PLOT_MARGINS['right'], top=segParam.PLOT_MARGINS['top'], bottom=segParam.PLOT_MARGINS['bottom'])

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show()
    modifierLabel = ""
    if normalize:
        modifierLabel += "norm_"
    if segParam.SKIP_OUTLIER_RUNS:
        modifierLabel += "skip_"
    fig.savefig(f"{saveFolder}Analysis/TimeSeries/CoM_Sq_{segParam.COM_MODE}_{modifierLabel}stdDev_initComp.png")
    plt.close(fig)
    print(f"Successfully plotted {modifierLabel}standard deviation of squared delta CoM {segParam.COM_MODE} for multiple architectures: {arcList} and initializations: {initList}")

def PlotSegregationComparison(runIndex: int) -> None:
    """Plots the segregation CoM trajectories of different architectures in a subplots figure.
    The architectures are specified in the AOI_COMPARE variable in Segregation_Parameters.py"""
    global architecture
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/big_bold.mplstyle")
    fig, ax = plt.subplots(nrows = segParam.NROWS, ncols = segParam.NCOLS, sharex = True, sharey = True)
    for n in range(len(segParam.AOI_COMPARE)):
        architecture = segParam.AOI_COMPARE[n]
        timeSteps, z_com_list = ReadData(runIndex)
        firstPassageTime, segTime = FindImprovedSegregationTime(z_com_list)
        # Coverting timesteps to units of TAU_0
        timeSteps, segTime, axisLabel = ScaleSegregationTimes(timeSteps, segTime)
        i = n % segParam.NROWS # row index
        j = n // segParam.NROWS # column index
        if segParam.NROWS == 1:
            axes_ij = ax[j]
        elif segParam.NCOLS == 1:
            axes_ij = ax[i]
        else:
            axes_ij = ax[i][j]
        PlotCoMDistance(axes_ij, timeSteps, z_com_list, rf"{axisLabel}$\tau_0$")
        if segTime != -1:
            PlotVerticalLine(axes_ij, segTime)
        axes_ij.set_title(architecture)
        axes_ij.set_xlabel("")
        axes_ij.set_ylabel("")
    ax[0].legend()
    fig.suptitle(f"Segregation Trajectory Comparison, {numberOfMonomers} monomers")
    fig.supxlabel(rf"Time Steps ({axisLabel}$\tau_0$)")
    yLabelModifier = ""
    if segParam.RESCALE_LENGTHS:
        yLabelModifier = r"/ $L$"
    fig.supylabel(r"$\Delta z_{CoM}$ %s" % (yLabelModifier))

    plt.show()
    folder = sysPaths.GetFolder(sysPaths.SEGREGATION, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)
    fig.savefig(f"{folder}segTrajectoryComparison.png")
    plt.close(fig)

def PlotComparisonInSinglePlot(ax: mpl.axes.Axes, runIndex: int, AOI_COMPARE: list[str] = segParam.AOI_COMPARE) -> None:
    """
    Plots the segregation time trajectories for architectures in the passed axis.
    Uses the default axis labels as used in PlotCoMDistance().
    """
    global architecture
    for arc in AOI_COMPARE:
        architecture = arc
        timeSteps, z_com_list = ReadData(runIndex)
        firstPassageTime, segTime = FindImprovedSegregationTime(z_com_list)
        timeSteps, segTime, axisLabel = ScaleSegregationTimes(timeSteps, segTime)
        PlotCoMDistance(ax, timeSteps, z_com_list, rf"{axisLabel}$\tau_0$", segParam.GetAliasArchitecture(architecture))
    # ax.set_xlabel("")
    # ax.set_ylabel("")
    # ax.set_xlim(left = -0.01, right = 0.4)
    # ax.legend()

def PlotOnlySinglePlotComparison(runIndex: int, showPlot: bool = False) -> None:
    """
    Plots the segregation trajectories for architectures indicated in AOI_COMPARE in a single plot for the passed run.
    """
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/big_bold.mplstyle")
    saveFolder = f"{sysPaths.GetFolder(sysPaths.SEGREGATION, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)}Segregation_Trajectory_Comparison/"
    # Making directory if it does not exist:
    pathlib.Path(saveFolder).mkdir(parents = True, exist_ok = True)
    fig, ax = plt.subplots()
    PlotComparisonInSinglePlot(ax, runIndex)
    # Setting Axis Labels
    ax.legend()
    ax.set_xlim(left = -0.1, right = 0.4)

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)

    # Saving figure
    fig.savefig(f"{saveFolder}seg_comparison_r{runIndex}{segParam.FIG_EXT}")
    plt.close(fig)
    print(f"Plotted Segregation Trajectory Comparison for {sysPaths.SPECIAL_SIMULATION} simulation Run {runIndex} for all architectures.")

def PlotOnlySinglePlotComparisonForAllRuns() -> None:
    """
    Plots the segregation trajectories for architectures indicated in AOI_COMPARE in a single plot for all runs.
    """
    for run in range(1, numberOfRuns + 1):
        PlotOnlySinglePlotComparison(run)
    
# Region wise CoM functions:
def ReadRegionData(folder: str, runIndex: int) -> Tuple[ArrayLike, List[ArrayLike]]: # type: ignore
    """Reads the CoM time series data for all regions and returns a tuple of the timesteps and the list of z coords of the CoM"""
    z_CoM_List = []
    timeSteps = []
    areTimeStepsSet = False #  A flag to indicate whether the timesteps have been read and set
    # Assuming all the regions will have the same number CoM Data calculated for the same timesteps
    for regionIndex in range(1, reg.NUMBER_OF_REGIONS + 1):
        filePath = f"{folder}run{runIndex}/com_reg{regionIndex}.dat"
        df = pd.read_csv(filePath)

        if not areTimeStepsSet:
            timeSteps = np.array(df.iloc[: , 0]) # reading the 1st colum: timesteps
            areTimeStepsSet = True
        z_CoM_List.append(np.array(df.iloc[: , 3])) # Reading the 4th column: z coord of CoM

    return timeSteps, z_CoM_List

def ReadAndPlotRegionData(showPlot: bool = False) -> None:
    """Invokes the ReadRegionData function and plots the CoM Time series for all the regions"""
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/big_bold.mplstyle")
    folder = GetFolder(numberOfMonomers, architecture)
    fig, ax = plt.subplots()
    # Setting colour scheme for regions:
    ax.set_prop_cycle(color = [ "#1D8AFE", '#00008B', "#FF3B3B", '#8B0000']) # orange: '#FFA500', dark red: '#8B0000'
    linestyles = ['-', '-', ':', ':']
    linewidths = [2, 2, 4, 4]

    timeSteps, z_CoM_List = ReadRegionData(folder, runIndex)


    # Converting time to units of tau_0
    timeSteps = timeSteps / segParam.TAU_0
    if not segParam.USE_LOGLOG_PLOT and not segParam.USE_SEMILOG_PLOT:
        timeSteps, axisLabel = pt.ConvertToScientificNotation(timeSteps)
    else:
        axisLabel = ""
    # Rescaling CoM positions with box length:
    yLabelModifier = r"($\sigma$)"
    if segParam.RESCALE_LENGTHS:
        z_CoM_List = pt.RescaleArraysByLength(z_CoM_List, boxLength)
        yLabelModifier = r"/ $L$"
    
    # Plotting:
    if segParam.USE_LOGLOG_PLOT:
        plottingFunction = ax.loglog
    elif segParam.USE_SEMILOG_PLOT:
        plottingFunction = ax.semilogx
    else:
        plottingFunction = ax.plot
    for i in range(len(z_CoM_List)):
        plottingFunction(timeSteps, z_CoM_List[i], label = reg.REGION_LABELS[i], linestyle = linestyles[i], linewidth = linewidths[i])
    # Labelling:
    ax.set_xlabel(rf"Time ({axisLabel}$\tau_0$)")
    ax.set_ylabel(r"$z_{COM}$ %s" % (yLabelModifier))
    if segParam.SHOW_TITLE:
        ax.set_title(f"CoM Time Series of various regions\n {numberOfPolymers} {architecture} polymer(s), Run {runIndex}")
    else:
        ax.text(0.75, 0.85, "%s\n%s\n%s" % (segParam.GetAliasArchitecture(architecture), rf"$N={numberOfMonomers}$", segParam.ORIENTATION_LABEL), transform = ax.transAxes, verticalalignment='center', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1))
    mpl.rcParams["legend.fontsize"] = 20
    leg = ax.legend(loc = 'center right')
    # Increasing legend line thickness:
    for line in leg.get_lines():
        line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())

    # Increasing x limit to accomodate legend and title:
    ax.set_xlim(right = 10 ** 8)

    if segParam.USE_MARGINS:
        fig.subplots_adjust(left=segParam.PLOT_MARGINS['left'], right=segParam.PLOT_MARGINS['right'], top=segParam.PLOT_MARGINS['top'], bottom=segParam.PLOT_MARGINS['bottom'])

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)
    # Saving the figure
    if segParam.USE_LOGLOG_PLOT:
        plotLabel = "loglog_"
    elif segParam.USE_SEMILOG_PLOT:
        plotLabel = "semilog_"
    else:
        plotLabel = ""
    fig.savefig(f"{folder}Analysis/Region_CoM_TimeSeries/{plotLabel}regCoM_r{runIndex}{segParam.FIG_EXT}")
    plt.close(fig)
    print(f"Successfully plotted region wise CoM Time Series for {architecture} Run {runIndex}")

def PlotRegionPairDistance(regionPairList: list[tuple], showPlot: bool = False) -> None:
    """Reads the centres of mass time series for all regions and then plots the distance
    between pairs of them as a function of time.
    The regionPairList argument is a list of pairs of region indices between which the distance is to 
    be found."""
    folder = GetFolder(numberOfMonomers, architecture)
    fig, ax = plt.subplots()
    timeSteps, z_CoM_List = ReadRegionData(folder, runIndex)

    # Converting time to units of tau_0
    timeSteps = timeSteps / segParam.TAU_0
    timeSteps, axisLabel = pt.ConvertToScientificNotation(timeSteps)
    # Rescaling CoM positions with box length:
    yLabelModifier = ""
    if segParam.RESCALE_LENGTHS:
        z_CoM_List = pt.RescaleArraysByLength(z_CoM_List, boxLength)
        yLabelModifier = r"/ $L$"

    # Plotting pairwise distance:
    for pair in regionPairList:
        # Validating pair of indices:
        for pairIndex in range(2):
            if not pair[pairIndex] in range(0, reg.NUMBER_OF_REGIONS):
                print(f"ERROR: The region index {pair[pairIndex]} of the pair {pair} is out of bounds! Terminating.")
                sys.exit(1)
        # Calculating distance:
        regionPairDistance = np.abs(z_CoM_List[pair[0]] - z_CoM_List[pair[1]])

        # Plotting:
        label = f"({reg.REGION_LABELS[pair[0]]})-({reg.REGION_LABELS[pair[1]]})"
        ax.plot(timeSteps, regionPairDistance, label = label)
    # Labelling:
    ax.set_xlabel(rf"Time steps ({axisLabel}$\tau_0$)")
    ax.set_ylabel(r"Pairwise $\Delta z_{CoM}$ %s" % (yLabelModifier))
    ax.set_title(f"Pairwise CoM Time Series of various regions\n {numberOfPolymers} {architecture} polymer(s), Run {runIndex}")
    leg = ax.legend()
    # Increasing legend line thickness:
    for line in leg.get_lines():
        line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())

    if showPlot:
        plt.show(block = True)
    # Saving the figure
    fig.savefig(f"{folder}Analysis/Region_CoM_TimeSeries/regDeltaCoM_r{runIndex}.png")
    plt.close(fig)
    print(f"Successfully plotted region wise CoM Time Series for {architecture} Run {runIndex}")
    

def PlotRegionDataForAllRuns() -> None:
    """Invokes ReadAndPlotRegionData() to plot regional CoM Time series for all runs"""
    global runIndex
    for runIndex in range(1, numberOfRuns + 1):
        ReadAndPlotRegionData()
    print(f"Plotted the Regional CoM Time Series for the {numberOfRuns} runs of {architecture}.")

def PlotCOMDistanceDistributionOnAxis(ax: mpl.axes.Axes, numberOfBins: int, runIndex: int = -1, label: str = "", readDistribution: bool = False, writeDistribution: bool = False) -> None:
    """
    Plots the probability distribution of the CoM distance between the two polymers on the given axis.
    If a runIndex is passed, the distribution for that run is plotted; otherwise the mean distribution over all runs is plotted.
    Args:
        - ax: The matplotlib axis on which to plot the distribution.
        - numberOfBins: The number of bins to use for the histogram.
        - runIndex: The index of the run to plot. If -1, the mean distribution over all runs is plotted.
        - label: The label to use for the plot legend.
        - readDistribution: If True, reads the distribution from a file instead of calculating it (given the file exists).
        - writeDistribution: If True, writes the calculated distribution to a file.
    """
    cutoffLabel = ""
    if segParam.USE_DISTRIBUTION_CUTOFF:
        cutoff_index = int(segParam.COM_DISTRIBUTION_CUTOFF / dataInterval)
        cutoffLabel = "cutoff_"
    folder = GetFolder(numberOfMonomers, architecture)
    if not readDistribution:
        if runIndex == -1: # averaging over all runs
            allDistances = []
            for i in range(1, numberOfRuns + 1):
                timeSteps, z_com_list = ReadData(i)
                distances = np.abs(z_com_list[0] - z_com_list[1])
                if segParam.USE_DISTRIBUTION_CUTOFF:
                    distances = distances[:cutoff_index]
                allDistances.extend(distances)
            distance = np.array(allDistances)
        else:
            timeSteps, z_com_list = ReadData(runIndex)
            distance = np.abs(z_com_list[0] - z_com_list[1])
            if segParam.USE_DISTRIBUTION_CUTOFF:
                distance = distance[:cutoff_index]

        # Rescaling the CoM distances with box length
        if segParam.RESCALE_LENGTHS:
            distance = distance / boxLength
        distribution, binEdges = np.histogram(distance, numberOfBins, density = True)
        binCenters = plotMonomerDensity.ShiftBinEdges(binEdges[:-1])
    else: # reading distribution from a file
        if runIndex == -1:
            filePath = f"{folder}run1/{cutoffLabel}COM_distance_distribution_all.csv"
        else:
            filePath = f"{folder}run{runIndex}/{cutoffLabel}COM_distance_distribution_r{runIndex}.csv"
        print(f"Reading CoM Distance Distribution from {filePath}")
        try:
            df = pd.read_csv(filePath)
        except FileNotFoundError:
            print(f"ERROR: The file {filePath} does not exist! Cannot read distribution. Terminating.")
            sys.exit(1)
        binCenters = np.array(df.iloc[:, 0])
        distribution = np.array(df.iloc[:, 1])
        
        
    # Rescaling the CoM distances with box length
    xLabelModifier = ""
    if segParam.RESCALE_LENGTHS:
        xLabelModifier = r"/ $L$"

    # Plotting:
    ax.plot(binCenters, distribution, marker = '.', linestyle = '--', label = label)
    ax.set_xlabel(rf"$\Delta z_{{CoM}}$ %s" % (xLabelModifier))
    ax.set_ylabel("Probability Density")

    # Writing distribution to a file if required, and if it is not already written down:
    if writeDistribution and not readDistribution: 
        if runIndex == -1:
            filePath = f"{folder}run1/{cutoffLabel}COM_distance_distribution_all.csv"
        else:
            filePath = f"{folder}run{runIndex}/{cutoffLabel}COM_distance_distribution_r{runIndex}.csv"
        with open(filePath, "w") as file:
            file.write(f"Bin Center, Probability Density\n")
            for i in range(len(binCenters)):
                file.write(f"{binCenters[i]}, {distribution[i]}\n")
        print(f"Wrote CoM Distance Distribution to {filePath}")

def PlotCOMDistributionRunsOnAxis(ax: mpl.axes.Axes, runIndices: list[int], fileLabel: str, label: str = "", numberOfBins: int = 500, readDistribution: bool = False, writeDistribution: bool = False) -> None:
    """
    Plots the probability distribution of the CoM distance between the two polymers for different runs of the same architecture on the same axis.
    The runs are specified in the runIndices argument.
    Args:
        - ax: The matplotlib axis on which to plot the distribution.
        - runIndices: The list of run indices to plot.
        - fileLabel: The label to use for the saved file.
        - label: The label to use for the plot legend.
        - numberOfBins: The number of bins to use for the histogram.
        - readDistribution: If True, reads the distribution from a file instead of calculating it (given the file exists).
        - writeDistribution: If True, writes the calculated distribution to a file.
    """
    cutoffLabel = ""
    if segParam.USE_DISTRIBUTION_CUTOFF:
        cutoff_index = int(segParam.COM_DISTRIBUTION_CUTOFF / dataInterval)
        cutoffLabel = "cutoff_"
    folder = GetFolder(numberOfMonomers, architecture)
    if not readDistribution:
        # averaging over all runs
        allDistances = []
        for i in runIndices:
            timeSteps, z_com_list = ReadData(i)
            distances = np.abs(z_com_list[0] - z_com_list[1])
            if segParam.USE_DISTRIBUTION_CUTOFF:
                distances = distances[:cutoff_index]
            allDistances.extend(distances)
        distance = np.array(allDistances)
        # Rescaling the CoM distances with box length
        if segParam.RESCALE_LENGTHS:
            distance = distance / boxLength
        distribution, binEdges = np.histogram(distance, numberOfBins, density = True)
        binCenters = plotMonomerDensity.ShiftBinEdges(binEdges[:-1])
    else: # reading distribution from a file
        filePath = f"{folder}run1/{cutoffLabel}COM_distance_distribution_{fileLabel}.csv"
        print(f"Reading CoM Distance Distribution from {filePath}")
        try:
            df = pd.read_csv(filePath)
        except FileNotFoundError:
            print(f"ERROR: The file {filePath} does not exist! Cannot read distribution. Terminating.")
            sys.exit(1)
        binCenters = np.array(df.iloc[:, 0])
        distribution = np.array(df.iloc[:, 1])
    
    # Rescaling the CoM distances with box length
    xLabelModifier = ""
    if segParam.RESCALE_LENGTHS:
        xLabelModifier = r"/ $L$"

    # Plotting:
    ax.plot(binCenters, distribution, marker = '.', linestyle = '--', label = label)
    ax.set_xlabel(rf"$\Delta z_{{CoM}}$ %s" % (xLabelModifier))
    ax.set_ylabel("Probability Density")

    # Writing distribution to a file if required, and if it is not already written down:
    if writeDistribution and not readDistribution: 
        filePath = f"{folder}run1/{cutoffLabel}COM_distance_distribution_{fileLabel}.csv"
        with open(filePath, "w") as file:
            file.write(f"Bin Center, Probability Density\n")
            for i in range(len(binCenters)):
                file.write(f"{binCenters[i]}, {distribution[i]}\n")
        print(f"Wrote CoM Distance Distribution to {filePath}")

def PlotCOMDistanceDistribution(numberOfBins: int = 500, runIndex: int = -1, showPlot: bool = False) -> None:
    """
    Plots the probability distribution of the CoM distance between the two polymers.
    If a runIndex is passed, the distribution for that run is plotted; otherwise the mean distribution over all runs is plotted.
    """
    folder = GetFolder(numberOfMonomers, architecture)
    if runIndex == -1:
        runIndexLabel = f"Mean over {numberOfRuns} runs"
        saveRunLabel = 'all'
    else:
        runIndexLabel = f"Run {runIndex}"
        saveRunLabel = f"r{runIndex}"
    cutoffLabel = ""
    if segParam.USE_DISTRIBUTION_CUTOFF:
        cutoffLabel = f", Cutoff = {int(segParam.COM_DISTRIBUTION_CUTOFF / segParam.TAU_0) :.1E} $\\tau_0$"
    # Plotting:
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/bold.mplstyle")
    fig, ax = plt.subplots()
    PlotCOMDistanceDistributionOnAxis(ax, numberOfBins, runIndex, readDistribution = segParam.READ_COM_DISTANCE_DISTRIBUTION, writeDistribution = segParam.WRITE_COM_DISTANCE_DISTRIBUTION)
    ax.set_title(f"COM Distance Distribution for two {architecture} polymers\n{runIndexLabel}, {numberOfMonomers} monomers, {sysPaths.SPECIAL_SIMULATION}{cutoffLabel}")

    if showPlot:
        plt.show(block = True)

    # Saving figure:
    saveFolder = f"{folder}Analysis/Probability_Distribution/"
    Path(saveFolder).mkdir(parents = True, exist_ok = True) # creating the folder if it does not exist
    fig.savefig(f"{saveFolder}COM_distance_distribution_{saveRunLabel}.png")
    plt.close(fig)
    print(f"Plotted CoM Distance Distribution for {architecture} {runIndexLabel}")

def PlotDistributionComparison(arcList: list[str], runIndex: int = -1, numberOfBins: int = 500, showPlot: bool = False) -> None:
    """
    Plots the probability distribution of the CoM distance between the two polymers for different architectures on the same axis.
    The architectures are specified in the arcList argument.
    """
    global architecture
    folder = GetFolder(numberOfMonomers, architecture)
    if runIndex == -1:
        runIndexLabel = f"Mean over {numberOfRuns} runs"
        saveRunLabel = 'all'
    else:
        runIndexLabel = f"Run {runIndex}"
        saveRunLabel = f"r{runIndex}"
    cutoffLabel = ""
    if segParam.USE_DISTRIBUTION_CUTOFF:
        cutoffLabel = f", Cutoff = {int(segParam.COM_DISTRIBUTION_CUTOFF / segParam.TAU_0) :.1E} $\\tau_0$"
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/bold.mplstyle")
    fig, ax = plt.subplots()
    for arc in arcList:
        architecture = arc
        PlotCOMDistanceDistributionOnAxis(ax, numberOfBins, runIndex, label = arc, readDistribution = segParam.READ_COM_DISTANCE_DISTRIBUTION, writeDistribution = segParam.WRITE_COM_DISTANCE_DISTRIBUTION)
        print(f"Plotted distribution for {arc}.")
    ax.set_title(f"COM Distance Distribution Comparison for two polymers\n{numberOfMonomers} monomers, {sysPaths.SPECIAL_SIMULATION}{cutoffLabel}")
    ax.legend()

    if showPlot:
        plt.show(block = True)
    # Saving figure:
    saveFolder = f"{folder}Analysis/Probability_Distribution/"
    Path(saveFolder).mkdir(parents = True, exist_ok = True) # creating the folder if it does not exist
    fig.savefig(f"{saveFolder}COM_distance_distribution_comp_{saveRunLabel}.png")
    plt.close(fig)
    print(f"Plotted CoM Distance Distribution Comparison for {arcList}, {numberOfMonomers} monomers")

def PlotDistributionRunComparison(numberOfBins: int = 500, showPlot: bool = False) -> None:
    """
    Plots the probability distribution of the CoM distance between the two polymers for different runs of the same architecture on the same axis.
    The runs are specified in the Segregation_Parameters.py script.
    """
    folder = GetFolder(numberOfMonomers, architecture)
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/bold.mplstyle")
    fig, ax = plt.subplots()
    segParam.UpdateGroupIndices(set(range(1, numberOfRuns + 1))) # updating the group indices to include all runs
    for n in range(len(segParam.GROUP_INDICES)):
        PlotCOMDistributionRunsOnAxis(ax, list(segParam.GROUP_INDICES[n]), fileLabel = f"grp{n+1}", label = segParam.GROUP_LABELS[n], numberOfBins = numberOfBins, readDistribution = segParam.READ_COM_DISTANCE_DISTRIBUTION, writeDistribution = segParam.WRITE_COM_DISTANCE_DISTRIBUTION)

    cutoffLabel = ""
    if segParam.USE_DISTRIBUTION_CUTOFF:
        cutoffLabel = f", Cutoff = {int(segParam.COM_DISTRIBUTION_CUTOFF / segParam.TAU_0) :.1E} $\\tau_0$"
    ax.set_title(f"COM Distance Distribution Comparison for two {architecture} polymers\n{numberOfMonomers} monomers{cutoffLabel}")
    ax.legend()

    if showPlot:
        plt.show(block = True)
    # Saving figure:
    saveFolder = f"{folder}Analysis/Probability_Distribution/"
    Path(saveFolder).mkdir(parents = True, exist_ok = True) # creating the folder if it does not exist
    fig.savefig(f"{saveFolder}COM_distance_distribution_run_grp_comp.png")
    plt.close(fig)
    print(f"Plotted CoM Distance Distribution Comparison for group of runs, {numberOfMonomers} monomers each.")

def TestOutliers():
    runIndices, segTimes = segParam.ReadSegregationTimes(numberOfMonomers, architecture, segParam.CRITERION_STRING)
    outlierIndices = segParam.FindOutlierIndices(segTimes)
    print(f"Outliers information for segregation times of {architecture} with {numberOfMonomers} monomers:")
    for i in outlierIndices:
        print(f" - Run {runIndices[i]}: Segregation Time = {segTimes[i]}")



# Performing operations:
if __name__ == "__main__":
    SetArchitectureAndMonomers()
    SetConstants()
    # print(f"Length of Box: {boxLength}")

    runIndex = GetOptionalArguments()
    if chosenMode == "timeseries":
        if runIndex == -1: # no run argument passed
        # Plotting comparison:
            if segParam.ONLY_TRAJECTORY_COMPARISON:
                PlotOnlySinglePlotComparisonForAllRuns()
                sys.exit(0)
            SegregateForAllRuns()
            if reg.USE_REGIONS:
                PlotRegionDataForAllRuns()
        else:
        # Plotting comparison:
            if segParam.ONLY_TRAJECTORY_COMPARISON:
                PlotOnlySinglePlotComparison(runIndex, showPlot = True)
                sys.exit(0)
            if not reg.USE_REGIONS:
                timeSteps, z_com_list = ReadData(runIndex)
                firstPassageTime, segTime = FindImprovedSegregationTime(z_com_list)
                PlotData(timeSteps, z_com_list, firstPassageTime, segTime, True)
                if(segTime == -1):
                    print(f"Run {runIndex} did not show segregation\n")
                else:
                    print(f"Segregation successful for run {runIndex}. Segregation Time = {segTime}, First Passage Time = {firstPassageTime}")
            else:
                ReadAndPlotRegionData(True)
                # PlotRegionPairDistance([(0, 2), (1, 3)], True)
            # PlotSegregationComparison(runIndex)
    elif chosenMode == "squared_timeseries":
        if segParam.PLOT_SIZE_COMPARISON:
            PlotSqDeltaCOMSizeComparison(segParam.SIZE_LIST, segParam.AOI_COMPARE, runIndex, True)
        elif segParam.PLOT_COMPARISON:
            PlotSqDeltaCOMArcComparison(segParam.AOI_COMPARE, True, runIndex)
        else:
            PlotSqDeltaCOMDisplacement(runIndex, True)
    elif chosenMode == "runs_squared_timeseries":
        if segParam.USE_GROUPS:
            PlotSqDeltaCOMGroupsComparison(True)
        else:
            PlotSquaredCOMRunComparison(list(range(1, 11)), True) # First ten independent runs
    elif chosenMode == "squared_stddev":
        if segParam.PLOT_INIT_COMPARISON:
            PlotSqCOMStdDevInitComparison(segParam.SPECIAL_SIMULATIONS, segParam.AOI_COMPARE, normalize=segParam.NORMALIZE_STD_DEV, showPlot = True)
        elif segParam.PLOT_COMPARISON:
            PlotSqCOMStdDevArcComparison(segParam.AOI_COMPARE, normalize = segParam.NORMALIZE_STD_DEV, showPlot = True)
        else:
            PlotSqCOMStdDev(normalize = segParam.NORMALIZE_STD_DEV, showPlot = True)
    elif chosenMode == "distribution":
        if segParam.PLOT_COMPARISON:
            if segParam.USE_GROUPS:
                PlotDistributionRunComparison(100, True)
            else:
                PlotDistributionComparison(segParam.AOI_COMPARE, runIndex, 100, True)
        else:
            PlotCOMDistanceDistribution(100, runIndex, True)
    elif chosenMode == "test_outliers":
        TestOutliers()
    print(f"Architecture: {architecture}")


# plotting CoM Distance:
# timeSteps, z_com_list = ReadData(runIndex)
# segTime = FindImprovedsegregationTIme(z_com_list)
# PlotCoMDistance(timeSteps, z_com_list, segTime)
