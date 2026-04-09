import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal
import sys
from typing import List, Tuple
import numpy.typing as npt
# importing the scripts that contains all the local paths
sys.path.append(f"../../Global_Scripts/System_File_Paths/")
import system_file_paths as sysPaths 
import Segregation_Parameters as segParam

# importing PlottingTools:
sys.path.append(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/")
import PlottingTools as pt

# polymer info:
architecture = "Arc2"
numberOfPolymers = 2
numberOfMonomers = 200

# plotting info:
plottingStyles = {"distribution": False, "scatter": False, "both": False} # dictionary with string-boolean key value pairs; a list of flags
plottingStyle = "scatter" # default style

# file info:
numberOfRuns = 50 # Number of runs that showed successful segregation
directory = sysPaths.SEGREGATION
initializationProcedure = ""

def VerifyPlottingStyle(plottingStyle: str) -> None:
    """Verifies whether the plotting style passed is a valid option that exists in the plotting styles list.
    Exits if it is invalid."""

    if not plottingStyle in plottingStyles.keys():
        print(f"The enetered plotting style was not recognized. Please choose from the available styles and try again: {plottingStyles.keys()}")
        exit(1)

def SetConstants():
    """Reads the argument(s) passed while invoking the script from the command line and sets it to the architecture and numberOfMonomers"""
    global architecture
    global numberOfMonomers
    global plottingStyle
    global initializationProcedure

    numberOfMandatoryArguments = 3
    if len(sys.argv) < numberOfMandatoryArguments:
        print("Not enough arguments specified while invoking the script! Please pass the number of monomers, the architecture, and the initialization procedure in the following format:")
        print("python /<path>/plotCoMDistribution.py <noOfMonomers> <architecture> <initializationProcedure>")
        print("The keyword 'all' can be entered instead of an architecture to plot a comparison plot for multiple architectures")
        print(f"Optionally, a fourth argument can be passed for the plotting style. Available options: {plottingStyles.keys()}")
        quit() # terminating the script
    
    architecture = sys.argv[2]
    numberOfMonomers = sys.argv[1]
    initializationProcedure = sys.argv[3]
    # checking if the number of monomers is valid:
    if not numberOfMonomers.isdigit():
        print(f"The entered numberOfMonomers {numberOfMonomers} cannot be converted to a number! Please provide a valid number")
        quit()

    numberOfMonomers = int(numberOfMonomers)

    if len(sys.argv) == numberOfMandatoryArguments + 2: # optional argument passed
        plottingStyle = sys.argv[numberOfMandatoryArguments + 1]
        VerifyPlottingStyle(plottingStyle)
    SetPlottingStyleFlagValues(plottingStyle)

def PrintStatistics(times: npt.ArrayLike, architecture: str = "") -> None:
    """Prints useful statistics of the passed segregation times: lower quartile, median, upper quartile to the console."""
    lower_quartile = np.percentile(times, 25)
    median = np.percentile(times, 50)
    upper_quartile = np.percentile(times, 75)

    # Printing statistics along with (approximate run number):
    print(f"{architecture} Statistics:")
    print(f"Lower Quartile:\t {lower_quartile} Tau_0")
    print(f"Median:\t\t {median} Tau_0")
    print(f"Upper Quartile:\t {upper_quartile} Tau_0")

def PlotHistogram(times: dict[int, int], ax, label: str, segregation_criterion: str, putAxisLabels: bool = True) -> None:
    """Plots the histogram of the passed times on the passes axis object"""
    timesList = list(times.values())
    if putAxisLabels:
        PrintStatistics(timesList) # Printing the statistics to console
        timesList, axisLabel = pt.ConvertToScientificNotation(timesList)
    ax.hist(timesList, bins = 25, label = label)
    if putAxisLabels:
        ax.set_xlabel(rf"Segregation Time ({axisLabel}$\tau_0$)")
        ax.set_ylabel("Frequency")
        ax.set_title(f"Distribution of Segregation Times for {numberOfPolymers} {segParam.GetAliasArchitecture(architecture)} polymers\n{numberOfMonomers} each, {len(timesList)} segregating runs, {segregation_criterion}")
    if len(label) > 0:
        ax.legend()

def PlotScatterPlot(times: dict[int, int], ax, label: str, segregation_criterion: str, putAxisLabels: bool = True) -> None:
    """Plots the scatter plot of the passed times (in Tau_0 units) on the passed axis object
    If putAxisLabels is True, the times are converted to scientific notation and axis labels are set accordingly."""
    x = list(times.keys())
    timesList = np.array(list(times.values()))
    if putAxisLabels:
        PrintStatistics(timesList) # Printing the statistics to console
        timesList, axisLabel = pt.ConvertToScientificNotation(timesList)
    ax.scatter(x, timesList, label = label)
    if putAxisLabels:
        ax.set_xlabel("Run Index")
        ax.set_ylabel(rf"Seg. Time ({axisLabel}$\tau_0$)")
        ax.set_title(f"Scatter Plot")
    if len(label) > 0:
        ax.legend()

def PlotBoxPlot(times: dict[int, int], ax, segregation_criterion: str, putAxisLabels: bool = True) -> None:
    """Plots a single box plot of the passed times on the axis object"""
    timesList = np.array(list(times.values()))
    if putAxisLabels:
        timesList, axisLabel = pt.ConvertToScientificNotation(timesList)
    ax.boxplot(timesList)
    if putAxisLabels:
        ax.set_title(f"Box Plot")
        # ax.set_ylabel(rf"Seg. Time ({axisLabel}$\tau_0$)")
    ax.tick_params(axis = 'x', bottom = False, labelbottom = False)

def read_csv_as_dict(filePath: str, header = 'infer', sep: str = ',') -> dict:
    """Reads a two-column CSV file at the given file path and returns each row as a key-value pair in a dictionary"""
    df = pd.read_csv(filePath, header = header, sep = sep)
    if len(df.columns) != 2:
        raise Exception(f"The csv file {filePath} has more than or less than 2 columns. Two columns were expected to convert the data to a dictionary.")
    # else:
    return dict(df.values)


def ReadSegTimes(folder: str) -> Tuple[dict[int, int], str]:
    """Reads the segregation times file and returns a tuple containing: 
    (i) the segregation times as a dictinary <index, time>, and
    (ii) A string stating the segregation criterion with which the segregation times were calculated"""
    filePath = ""
    # if sysPaths.IsCylinderInfinite():
    #     segregation_criterion = segParam.INF_CRITERION_STRING
    # else:
    # Using the usual axis-length-based CoM difference to discern segregation
    segregation_criterion = segParam.CRITERION_STRING
    filePath = f"{folder}segregationTimes_{segregation_criterion}.csv"
    df = pd.read_csv(filePath)
    # Getting segregation string from the last column header:
    read_segregation_criterion = pt.ExtractSubstringInParenthesis(df.columns.values.tolist()[-1])
    if read_segregation_criterion != segregation_criterion:
        print("ERROR: The segregation criterion string in the file name and the file header are different!")
        print(f"File: {filePath}")
        sys.exit(1)

    return read_csv_as_dict(filePath), segregation_criterion

def GetStatisticsLabel(times: npt.ArrayLike) -> str:
    """Returns the mean, std, and relative error as a single string to add to the graph"""
    mean = np.mean(times)
    std = float(np.std(times))
    relativeError = std / mean # The fractional error for the segregation time
    label = f"Mean = {Decimal(mean):.2E}\n Std = {Decimal(std):.2E}\nRelative Error = {Decimal(relativeError):.2E}"
    return label

# TODO: update this list of boolean flags to a dictionary if more styles are introduced
def SetPlottingStyleFlagValues(plottingStyle: str) -> None:
    """Sets the boolean flags in the dictionary to indicate which plots to make."""
    if plottingStyle == "both":
        for key in plottingStyles.keys():
            plottingStyles[key] = True # Setting all flags to True
        plottingStyles["both"] = False # resetting the 'both' flag since it will not be functional
    elif plottingStyle in plottingStyles.keys():
        plottingStyles[plottingStyle] = True
    else:
        print(f"The plotting style {plottingStyle} is not recognized!")
        sys.exit(1)

def RescaleTimesDict(timesDict: dict[int, int], rescaleFactor: float) -> dict[int, int]:
    """Rescales the times of a times dictionary by dividing by the rescale factor."""
    newDict = {}
    for run in timesDict.keys(): # scaling by Tau_0
        newDict[run] = timesDict[run] / rescaleFactor # time in units of Tau_0
    return newDict
        

def ReadAndPlotSegTimes(folder: str, plottingStyle: str, showPlot: bool = False) -> None:
    """Reads the segregation times file and plots the data"""

    # Setting the plotting style:
    mpl.style.use(f'{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/bold.mplstyle')

    times, segregation_criterion = ReadSegTimes(folder)
    times = RescaleTimesDict(times, segParam.TAU_0)
    label = GetStatisticsLabel(list(times.values()))

    # Plotting Style:
    VerifyPlottingStyle(plottingStyle)
    SetPlottingStyleFlagValues(plottingStyle)

    # plotting:
    fileNames = {"scatter": "segTimes.png", "distribution": "segTimeDistribution.png"} # Dictionary to return the filename of the figure to be saved
    plottingFunction = {"scatter": PlotScatterPlot, "distribution": PlotHistogram} # Dictionary to return the delegate to be called to plot the times
    for key in plottingStyles.keys():
        if plottingStyles[key]:
            fig, ax = plt.subplots()
            plottingFunction[key](times, ax, label = label, segregation_criterion = segregation_criterion)
            fig.savefig(f"{folder}{fileNames[key]}")
            print(f"Succesfully plotted Segregation Time {key} plot for {architecture}")
    
    # plt.xticks(fontweight = 'bold')
    if showPlot:
        plt.show(block = True)

def ReadAndPlotSegTimesWithBoxPlot(folder: str, plottingStyle: str, showPlot: bool = False) -> None:
    """Plots the segregation times distribution and the box plot side by side for an architecture"""
    mpl.style.use(f'{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/subplots_big_bold.mplstyle')

    times, segregation_criterion = ReadSegTimes(folder)
    times = RescaleTimesDict(times, segParam.TAU_0)
    if segParam.SHOW_TITLE:
        label = GetStatisticsLabel(list(times.values()))
    else:
        label = ""

    # Plotting Style:
    VerifyPlottingStyle(plottingStyle)
    SetPlottingStyleFlagValues(plottingStyle)

    fig, axes = plt.subplots(nrows = 1, ncols = 2, sharey = True, gridspec_kw = {"width_ratios": [3, 1]})
    
    # Plotting distribution:
    if plottingStyle == "distribution":
        PlotHistogram(times, axes[0], label = label, segregation_criterion = segregation_criterion, putAxisLabels = True)
    elif plottingStyle == "scatter":
        PlotScatterPlot(times, axes[0], label = label, segregation_criterion = segregation_criterion, putAxisLabels = True)
    else:
        print("Style not found: '%s'. Please use a valid specific style." % (plottingStyle))
        exit(1)
    
    # Plotting the box plot:
    PlotBoxPlot(times, axes[1], True)
    
    if showPlot:
        plt.show(block = True)

    # Saving figure:
    fig.savefig(f"{folder}segTimesBoxPlot.png")

def GetAxis(axes, nrows: int, ncols: int, i: int, j: int) -> mpl.axes._axes.Axes:
    """Returns the ijth axis matrix element of the subplots axes object"""
    if nrows == 1:
        axes_ij = axes[j]
    elif ncols == 1:
        axes_ij = axes[i]
    else:
        axes_ij = axes[i][j]
    return axes_ij

def ReadAndPlotGroupedSegTimes(folder: str, plottingStyle: str, showPlot: bool = False) -> None:
    """Reads the segregation times from the folder and plots them in multiple groups, as defined in the segParam file"""
    # Setting the plotting style:
    mpl.style.use(f'{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/subplots_big_bold.mplstyle')

    # Reading times:
    rawTimes, segregation_criteria = ReadSegTimes(folder)
    rawTimes = RescaleTimesDict(rawTimes, segParam.TAU_0) # time in units of Tau_0
    
    timesList, axisLabel = pt.ConvertToScientificNotation(np.array(list(rawTimes.values())))
    times = dict(zip(rawTimes.keys(), timesList)) # new dictionary with the rescaled times

    # Plotting Style:
    VerifyPlottingStyle(plottingStyle)
    SetPlottingStyleFlagValues(plottingStyle)

    # Setting groups
    segParam.UpdateGroupIndices(set(range(1, numberOfRuns+1)))
    # debugging:
    # print(segParam.GROUP_INDICES)
    group_times = []
    labels = [] # a list of statistical labels for each group
    for n in range(segParam.NGROUPS): # Adding the group times from the indices in GROUP_INDICES
        # group_times.append([times[i-1] for i in segParam.GROUP_INDICES[n]])
        groupDict = {}
        rawTimesList = [] # A list for the raw (only sacled by Tau_0) seg times for each group
        for run in segParam.GROUP_INDICES[n]:
            if run in times.keys(): # ignoring runs that do not occur in csv file; unsegregated 
                rawTimesList.append(rawTimes[run])
                groupDict[run] = times[run] # adding (run index, time) to the group dictionary
        labels.append(GetStatisticsLabel(rawTimesList))
        group_times.append(groupDict) # Adding group dictionary for each group

    

    # Plotting:
    nrows = segParam.NROWS
    ncols = segParam.NCOLS
    fig, axes = plt.subplots(nrows = nrows, ncols = ncols, sharey = True)
    for n in range(segParam.NGROUPS):
        # Getting the ijth element of the subplot axes:
        i = n % nrows
        j = n // nrows
        axis = GetAxis(axes, nrows, ncols, i, j)
        # Plotting distribution:
        if plottingStyle == "distribution":
            PlotHistogram(group_times[n], axis, label = labels[n], segregation_criterion = segregation_criteria, putAxisLabels = False)
        elif plottingStyle == "scatter":
            PlotScatterPlot(group_times[n], axis, label = labels[n], segregation_criterion = segregation_criteria, putAxisLabels = False)
        else:
            print("Style not found: '%s'. Please use a valid specific style." % (plottingStyle))
            sys.exit(1)
        axis.set_title(f"{segParam.GROUP_LABELS[n]} ({numberOfMonomers})") # Setting title on each subplot
    # Setting sub titles, xlabels, ylabels:
    if segParam.SHOW_TITLE:
        fig.suptitle("Segregation Times for various groups of %s\n %i monomers, Seg. criterion: %s" % (architecture, numberOfMonomers, segregation_criteria))
    if plottingStyles["scatter"]:
        fig.supxlabel("Run Index")
        fig.supylabel(rf"Segregation Time ({axisLabel}$\tau_0$)")
    elif plottingStyles["distribution"]:
        fig.supxlabel(rf"Segregation Time ({axisLabel}$\tau_0$)")
        fig.supylabel("Frequency")
    
    if showPlot:
        plt.show(block = True)

    # Saving figure:
    fig.savefig(f"{folder}segTimeComparison{segParam.FIG_EXT}")
    plt.close(fig)

    # Printing grouped segregation times to console in csv format:
    for n in range(segParam.NGROUPS):
        print(f"Group {segParam.GROUP_LABELS[n]}")
        print("Run Index, Segregation Time")
        for index in group_times[n]:
            print(f"{index}, {rawTimes[index]} Tau_0")

def GetFolder(architecture: str) -> str:
    """Returns the path to the analysis folder where the segregation times are stored for a particular architecture"""
    prefix = sysPaths.GetFolder(directory, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)
    return f"{prefix}{architecture}/Analysis/"

def GetArchiveFolder(segregationCriterion: str, architecture: str) -> str:
    """Returns the path of the segregatio times that have been archived as a particular segregation criterion"""
    return f"{directory}Analysis/Segregation_Times_Images/b{numberOfMonomers}/{segregationCriterion}/{architecture}/"

def GetPreviousAttemptsFolder(architecture: str, simulationName: str) -> str:
    """Returns the path of the segregation times file in a particular backed-up smulation folder"""
    return f"{directory}b{numberOfMonomers}/Previous_Attempts/{simulationName}/{architecture}/Analysis/"

def PlotCollageWithStyle(architectures: list[str], folder: str, plottingStyle: str, showPlot: bool = False, nrows: int = 2, ncols: int = 3) -> None:
    """Plots the segregation times for all the different architectures passed as a collage"""

    # Setting the plotting style:
    mpl.style.use(f'{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/subplots_bold.mplstyle')
    
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, sharey = True, sharex = True)

    # reading segregation times:
    times_array = [] # A list to store just the arrays of segregation times for each architecture
    timesDict_array = [] # A list to store the times dictionaries for each architecture
    labels = []
    segregation_criterion = "" # The criterion with which all segregations were supposedly carried out
    for n in range(len(architectures)):
        times, criterion = ReadSegTimes(GetFolder(architectures[n]))
        # Verifying if all segregation criteria are the same:
        if len(segregation_criterion) == 0: # not initialized
            segregation_criterion = criterion
        else:
            if segregation_criterion != criterion:
                print(f"ERROR: All the segregation criteria for the architecures are not the same! {segregation_criterion} and {criterion} for {architectures[n]} cannot be compared.")
                sys.exit(1)

        times = RescaleTimesDict(times, segParam.TAU_0)
        label = GetStatisticsLabel(list(times.values()))
        times_array.append(np.array(list(times.values())))
        timesDict_array.append(times)
        labels.append(label)
    
    times_array, axisLabel = pt.ConvertMultipleArraysToScientificNotation(times_array)
    for i in range(len(timesDict_array)):
        timesDict_array[i] = dict(zip(timesDict_array[i].keys(), times_array[i])) # Creating new dicts with rescaled times

    aliasArchitectures = segParam.GetAliasArchitectureList(architectures)
    for n in range(len(architectures)):
        i = n % nrows # row index
        j = n // nrows # column index
        if nrows == 1:
            axes_ij = axes[j]
        elif ncols == 1:
            axes_ij = axes[i]
        else:
            axes_ij = axes[i][j]
        if plottingStyle == "distribution":
            PlotHistogram(timesDict_array[n], axes_ij, label = labels[n], segregation_criterion = segregation_criterion, putAxisLabels = False)
        elif plottingStyle == "scatter":
            PlotScatterPlot(timesDict_array[n], axes_ij, label = labels[n], segregation_criterion = segregation_criterion, putAxisLabels = False)
        else:
            print("Style not found: '%s'. Please use a valid specific style." % (plottingStyle))
            exit(1)
        axes_ij.set_title("%s" % (aliasArchitectures[n]))
        # axes[i][j].tick_params(axis = 'both', labelleft = True, labelbottom = True)
        # zooming in:
        # axes[i][j].set_ylim(0, 5*10**6)
    # fig.tight_layout()
    # fig.subplots_adjust(hspace = 0.4)
    
    fig.suptitle("Segregation Times for various architectures, %i monomers, Seg. criterion: %s" % (numberOfMonomers, segregation_criterion))
    if plottingStyles["scatter"]:
        fig.supxlabel("Run Index")
        fig.supylabel(rf"Segregation Time ({axisLabel}$\tau_0$)")
    elif plottingStyles["distribution"]:
        fig.supxlabel(rf"Segregation Time ({axisLabel}$\tau_0$)")
        fig.supylabel("Frequency")
    
    if showPlot:
        plt.show(block = True)
    fig.savefig(f"{sysPaths.GetFolder(sysPaths.SEGREGATION, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)}Segregation_Times_subset.png")
    # fig.savefig(f"{sysPaths.POLYMER_PHYSICS}/LAMMPS_runs/cluster-data/new_segregation/Analysis/Segregation_Times_subset.png")

# TODO: Use architecture aliases
if __name__ == "__main__":
    # Setting architecture:
    SetConstants()
    if architecture != 'all':
        # Plotting:
        folder = GetFolder(architecture)
        ReadAndPlotSegTimes(folder, plottingStyle, True)
        # ReadAndPlotGroupedSegTimes(folder, plottingStyle, True)
        # ReadAndPlotSegTimesWithBoxPlot(folder, plottingStyle, True)
    else:
        # Plotting collage
        architectures = segParam.AOI_LIST
        folder = f"{directory}Analysis/"
        PlotCollageWithStyle(architectures, folder, plottingStyle, showPlot= True, nrows = segParam.NROWS, ncols = segParam.NCOLS)
