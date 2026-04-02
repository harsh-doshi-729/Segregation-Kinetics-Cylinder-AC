# This script is used to plot the distance between two loci (monomers) in a polymer system 
# as a function of time

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

# polymer information:
numberOfPolymers = 2
numberOfMonomers = 200 # each 
architecture = ""
numberOfRuns = segParam.NUMBER_OF_RUNS
runIndex = -1
indexPair = [] # A pair of the monomer indices representing the loci between which the distance is to be plotted

# Argument information:
numberOfAdditionalArguments ={"new_segregation": 1, "create_initial_states": 1, "custom": 1, "comparison": 0, "statistics": 0}
keywordInstructions = ["1 additional argument; an integer argument for the run index",
                        "1 additional argument; an integer argument for the run index",
                        "1 additional argument; a string argument for the custom directory path where the data file resides",
                        "0 additional arguments;",
                        "1 additional arguments; a string argument for the custom directory path where all the run folders reside"]
chosenKeyword = ""
directoryPath = ""
columnToBePlotted = ""

def PrintArgumentInstructions() -> None:
    """
    Prints the instructions for the directory keyword arguments required to run this script.
    """
    print("The other arguments depend on the directory keyword:")
    for keyword, instruction in zip(numberOfAdditionalArguments.keys(), keywordInstructions):
        print(f"For '{keyword}': {instruction}")

def ParseAdditionalArguments(numberOfMandatoryArguments: int) -> None:
    """
    Parses the arguments after the directory keyword: including the keyword arguments and the optional arguments.
    """
    global runIndex
    global directoryPath
    global columnToBePlotted

    if chosenKeyword == "new_segregation" or chosenKeyword == "create_initial_states":
        try:
            runIndex = int(sys.argv[numberOfMandatoryArguments + 1])
        except ValueError:
            print(f"ERROR: The additional argument {sys.argv[numberOfMandatoryArguments + 1]} could not be interpreted as an integer for the run index.")
            sys.exit(1)
            customFolderPath = sys.argv[numberOfMandatoryArguments + 1]
        if chosenKeyword == "new_segregation":
            baseFolder = sysPaths.NEW_SEGREGATION
        else:
            baseFolder = sysPaths.CREATE_INITIAL_STATES 
        directoryPath = f"{sysPaths.GetFolder(baseFolder, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)}{architecture}/run{runIndex}/"
    elif chosenKeyword == "custom" or chosenKeyword == "statistics":
        customFolderPath = sys.argv[numberOfMandatoryArguments + 1]
        directoryPath = customFolderPath

    # Optional arguments, if any are passed:
    if len(sys.argv) > numberOfMandatoryArguments + numberOfAdditionalArguments[chosenKeyword] + 1:
        columnToBePlotted = sys.argv[numberOfMandatoryArguments + numberOfAdditionalArguments[chosenKeyword] + 1]
    else:
        columnToBePlotted = "distance" # default value for the column 

def SetConstants() -> None:
    """Reads the argument(s) while invoking the script and sets the values of the global variables."""
    global numberOfMonomers
    global architecture
    global chosenKeyword
    global indexPair
    global directoryPath

    numberOfMandatoryArguments = 5
    if len(sys.argv) < numberOfMandatoryArguments + 1:
        print("Not enough arguments specified while invoking the script! Please pass the name of the number of monomers, architecture, the two indices of the monomers to plot the distance between, the directory keyword, and any additional required arguments in the following format:")
        print("python /<path>/script.py <noOfMonomers> <architecture> <monomer-index1> <monomer-index2> <directory-keyword> <other-args>")
        print("The directory keyword can be one of the following: " + ", ".join(numberOfAdditionalArguments.keys()))
        PrintArgumentInstructions()
        print("An optional argument for the column name indicated the column to be plotted can be passed the above arguments.")
        sys.exit(1)
    else:
        try:
            numberOfMonomers = int(sys.argv[1])
        except ValueError:
            print(f"ERROR: The argument {sys.argv[1]} could not be interpreted as a number of monomers!")
            sys.exit(1)
        
        architecture = sys.argv[2]

        try:
            indexPair = [int(sys.argv[3]), int(sys.argv[4])]
        except ValueError:
            print(f"ERROR: The arguments {sys.argv[3]} and {sys.argv[4]} could not be interpreted as indices of monomers!")
            sys.exit(1)

        chosenKeyword = sys.argv[5]
        if chosenKeyword not in numberOfAdditionalArguments.keys():
            print(f"ERROR: The directory keyword '{chosenKeyword}' is not recognized! Please choose from {numberOfAdditionalArguments.keys()}.")
            PrintArgumentInstructions()
            sys.exit(1)

        # Reading additional arguments based on the directory keyword:
        if len(sys.argv) <= numberOfMandatoryArguments + numberOfAdditionalArguments[chosenKeyword]:
            print(f"ERROR: Not enough additional arguments provided for the '{chosenKeyword}' directory keyword! Expected {numberOfAdditionalArguments[chosenKeyword]} additional arguments.")
            sys.exit(1)

        ParseAdditionalArguments(numberOfMandatoryArguments)

def GetReadFileName(indexPair: list[int]) -> str:
    """
    Returns the name of the distance data file based on the monomder indices chosen.
    """
    if len(indexPair) != 2:
        print(f"ERROR: While getting distance file name, the list of integers must have only two entries. Found {len(indexPair)}.")
        sys.exit(1)
    return f"distance-{indexPair[0]}-{indexPair[1]}.csv"

def PlotDistanceData(ax: mpl.axes.Axes, timeSteps: ArrayLike, y_quantity: ArrayLike, label: str = "") -> None:
    """
    Plots the passed distance data againts the passed timeSteps in the passed axis.
    Scales the timesteps by Tau_0 and converts it to scientific notation. Sets the x-axis label.
    """
    # Scaling timesteps:
    timeSteps = timeSteps / segParam.TAU_0
    timeSteps, axisLabel = pt.ConvertToScientificNotation(timeSteps)

    ax.plot(timeSteps, y_quantity, label = label)
    ax.set_xlabel(rf"Time Steps ({axisLabel}$\tau_0$)")

def ReadAndPlotDistanceData(showPlot: bool = False) -> None:
    """Reads the distance data from the specified file and plots it."""
    global columnToBePlotted
    # Reading the distance data:
    filePath = f"{directoryPath}{GetReadFileName(indexPair)}"
    try:
        df = pd.read_csv(filePath, skipfooter = 2, engine='python') # Skipping the last two statistics lines
    except FileNotFoundError:
        print(f"ERROR: The file {filePath} does not exist! Please check the path and try again.")
        sys.exit(1)

    timeSteps = np.array(df.iloc[:, 0])
    try:
        quantity = np.array(df[f" {columnToBePlotted}"])
    except KeyError:
        print(f"ERROR: The column ' {columnToBePlotted}' could not be found in the distance data file.")
        sys.exit(1)


    # Plotting the distance data:
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/bold.mplstyle")
    fig, ax = plt.subplots()
    PlotDistanceData(ax, timeSteps, quantity, label = columnToBePlotted)
    ax.set_ylabel("Distance")
    if chosenKeyword == "custom":
        secTitle = "Custom File Path"
    else:
        secTitle = f"Architecture: {architecture}, Run Index: {runIndex}"
    ax.set_title(f"Distance Time Series for Monomers {indexPair[0]} and {indexPair[1]}\n{secTitle}")
    plt.legend()
    
    if showPlot:
        plt.show(block = True)

    # Save the plot
    outputFileName = f"{directoryPath}{columnToBePlotted}_time_series_{indexPair[0]}_{indexPair[1]}{segParam.FIG_EXT}"
    fig.savefig(outputFileName)
    plt.close(fig)
    print(f"Distance time series plotted for {secTitle}.")

def PlotComparison() -> None:
    """
    Reads multiple distance data files (paths taken for Segregation_Parameters.py) and plots their data in a single figure.
    """
    fileName = GetReadFileName(indexPair)
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/bold.mplstyle")
    fig, ax = plt.subplots()
    for key in segParam.COMPARE_CUSTOM_PATHS.keys():
        # Reading the data file:
        filePath = f"{segParam.COMPARE_CUSTOM_PATHS[key]}{fileName}"
        df = pd.read_csv(filePath, skipfooter = 2, engine = 'python')
        timeSteps = np.array(df.iloc[:, 0])
        try:
            y_quantity = np.array(df[f" {columnToBePlotted}"])
        except KeyError:
            print(f"ERROR: The column ' {columnToBePlotted}' could not be found in the distance data file.")
            sys.exit(1)
        PlotDistanceData(ax, timeSteps, y_quantity, label = key)
    # Setting rest of axis elements and labels
    ax.set_ylabel("Distance")
    ax.set_title(f"Comparison of {columnToBePlotted} time-series")
    plt.legend()

    plt.show(block = True)
    plt.close(fig)

def GetStatistics(filePath: str) -> tuple[list[str], list[float], list[float]]:
    """
    Reads the last few lines for a distance data file and extracts the statistics of the time series for each distance calculated.
    Args:
        filePath: The path to the z_extent.txt file
    Returns:
        A tuple of a string list and two float lists. The string list contains the names of columns for which the statistics are calculated.
        The two lists are the means and standard deviations of the distances in each column.
    """
    # Reading the header / column names:
    with open(filePath, 'r') as file:
        columnNames = file.readline().strip().split(',')
        columnNames = [name.strip() for name in columnNames]  # Stripping whitespace from column names
        columnNames = columnNames[1:] # Skipping the first column (timeSteps)

    # Reading the last few statistics lines
    numberOfStatisticsLines = 2 # A header info line + a line for each region
    lastLines = pt.ReadLastLines(filePath, numberOfStatisticsLines)
    # Debugging:
    # print(f"Last lines (len = {len(lastLines)}) of z_extent:\nf{lastLines}")
    means = []
    stds = []
    try:
        meanString = lastLines[0].split(':')[-1] # Extracting the string after the colon
        words = meanString.split(',')
        means = [float(word.strip()) for word in words] # Converting the words to floats
        stdString = lastLines[1].split(':')[-1] # Extracting the string after the colon
        words = stdString.split(',')
        stds = [float(word.strip()) for word in words] # Converting the words to floats

    except ValueError:
        print(f"Some of the words {words} could not be converted to a float value.")
        sys.exit(1)

    return columnNames, means, stds

def CalculateDistanceSampleStatistics(runFolder: str) -> None:
    """
    For each of the z_extent.txt files for different runs, this function calculates mean of the sample and the standard error (standard deviation of the means) for each region.
    Args:
        runFolder: the folder where the run folders are present; each run folder should have a z_extent.txt file.
    """
    regionMeansList = []
    regionSTDsList = []
    fileName = GetReadFileName(indexPair)
    for run in range(1, segParam.NUMBER_OF_RUNS + 1):
        filePath = f"{runFolder}run{run}/{fileName}"
        columnNames, means, stds = GetStatistics(filePath)
        regionMeansList.append(means)
        regionSTDsList.append(stds)
    meanDF = pd.DataFrame(regionMeansList)
    stdDF = pd.DataFrame(regionMeansList)
    # Printing statistics:
    print(f"Statistics for distances between {indexPair[0]} and {indexPair[1]} monomers over {segParam.NUMBER_OF_RUNS} runs in {runFolder}:")
    for i in range(len(columnNames)):
        regionMeans = meanDF.iloc[:, i]
        standardError = np.std(regionMeans)
        stdErrSTD = standardError / 2 / (segParam.NUMBER_OF_RUNS - 1) # According to the scaled chi-squared distribution, apparently
        print(f"{columnNames[i]}: Mean = {np.mean(regionMeans):.6f}, SD of Means = {standardError:.6f}, SE of SD = {stdErrSTD:.6f}")

if __name__ == "__main__":
    SetConstants()
    if chosenKeyword == "comparison":
        PlotComparison()
    elif chosenKeyword == "statistics":
        CalculateDistanceSampleStatistics(directoryPath)
    else:
        ReadAndPlotDistanceData(True)
