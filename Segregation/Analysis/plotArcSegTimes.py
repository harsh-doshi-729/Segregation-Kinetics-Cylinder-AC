# This script plots the summary of distributions of segregation times across various architectures in one plot
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

# modules I have written:
sys.path.append(f"{sysPaths.CREATE_INITIAL_STATES}Scripts/")
import plotAverageCoMDistance
import plotSegTimeDistribution

# Global variables
MARGIN = 0.1 # The margin between the last annoated flier and the edge of the plot area
numberOfMonomers = 200
special_simulation = "" # The name of the particular algorithm used to generate mixed states and segregate them
showFliersInformation = True # A flag to indicate whether the box plots should have annotated text to show number of outliers/fliers

def GetAcceptedSpecialSimulationNames(numberOfMonomers: int) -> list[str]:
    """Returns a list of accepted special simulation names that can be used to make the plot"""
    filePath = f"{sysPaths.NEW_SEGREGATION}b{numberOfMonomers}/Previous_Attempts/" # The path to the Previous_Attempts folder 
    # which stores all the special simulation folders
    previousAttemptsPath = Path(filePath)

    # Checking existence of Previous_Attempts:
    if not previousAttemptsPath.exists:
        print(f"The directory {filePath} does not exist! This is where the special simulations should be. Terminating.")
        sys.exit(1)
    
    simulations = [dir.name for dir in previousAttemptsPath.iterdir() if dir.is_dir()]
    return simulations

def SetConstants():
    """Reads the argument(s) passed while invoking this script and stores them in the global variables"""
    global numberOfMonomers
    global special_simulation

    numberOfMandatoryArguments = 2
    if len(sys.argv) < numberOfMandatoryArguments + 1:
        print(f"Not enough arguments passed! Please pass the number of monomers and the special simulation name while invoking the script in the format specified below:")
        print("python plotArcSegTimes.py <noOfMonomers> <special_simulation>")
        print(f"Accepted simulation names: {GetAcceptedSpecialSimulationNames(numberOfMonomers)}")
        quit()
    
    try:
        numberOfMonomers = int(sys.argv[1])
    except ValueError:
        print(f"The passed argument {sys.argv[1]} could not be converted to an integer. Please enter a valid number.")
        quit()
    
    special_simulation = sys.argv[2]
    # Verifying if the special_simulation is in the accepted names:
    acceptedNames = GetAcceptedSpecialSimulationNames(numberOfMonomers)
    if not special_simulation in acceptedNames:
        print(f"The simulation name {special_simulation} could not be recognized. Please enter one from the following list:")
        print(f"{acceptedNames}")
        sys.exit(1)
    


def IsSubdirectoryAnArchitecture(subdirectoryName):
    """Returns True if the passed subdirectory is named as an architecture (starts with 'Arc'), otherwise returns false"""
    # checking is the name starts with 'Arc':
    if subdirectoryName[:3] == "Arc":
        return True
    else: 
        return False

def GetArchitectureList(directory: str):
    """Reads the names of the subdirectories in the passed directory that correspond to architecture names. Returns the list of names"""
    path = Path(directory)
    if not path.exists:
        print(f"The directory {directory} does not exist!")
        quit(1)
    subdirectories = path.iterdir()

    architectures = []
    for subdirectory in subdirectories:
        name = subdirectory.name
        if subdirectory.is_dir() and IsSubdirectoryAnArchitecture(name):
            architectures.append(name)
    
    return architectures

def GetSegTimeDistribution(folder: str, architecture: str) -> tuple[list[int], str]:
    """Reads the segregation times for  an architecture and returns a tuple of two items:
    (i) The list of segregation times for the different runs of the architecture specified in the folder destination
    (ii) The segregation criterion string that was used to find segregation times."""
    folderPath = f"{folder}{architecture}/Analysis/"
    timesDict, segregation_criterion = plotSegTimeDistribution.ReadSegTimes(folderPath)
    segTimes = list(timesDict.values())
    return segTimes, segregation_criterion

def PlotMeans(architectures, folder: str, normalize: bool):
    """Plots the means of the segregation times of different architectures passed in the list and saves the image at the location of the folder.
    normalize: if True, the segregation time means are 'normalized' by the average distance the CoMs must move further to be called segregated.
    This normalization is done with a quantity that signifies how easy it is for two polymers to segregate from a given mixed state separation."""

    if normalize:
        plotAverageCoMDistance.SetConstants() # initializing the variables in this module
        distancesToSegregate = []
    means = []
    stds = []
    segregation_criterion = "" # The criterion with which all segregations were supposedly carried out
    for arc in architectures:
        segTimes, criterion = GetSegTimeDistribution(folder, arc)
        if len(segregation_criterion) == 0: # not initialized
            segregation_criterion = criterion # initializing with criterion of first architecture
        else:
            if segregation_criterion != criterion: # comparing criteria of different architectures
                print(f"ERROR: All the segregation criteria for the architecures are not the same! {segregation_criterion} and {criterion} for {arc} cannot be compared.")
                sys.exit(1)
        dist = np.array(segTimes)
        mean = np.mean(dist)
        std = np.std(dist)

        if normalize:
            distanceToSegregate = plotAverageCoMDistance.GetDistanceToSegregation(arc) # this parameter is lower for mixed state that have a higher average CoM distance
            distancesToSegregate.append(distanceToSegregate)
            mean = mean / distanceToSegregate
            std = std / distanceToSegregate

        means.append(mean)
        stds.append(std)

    fig, ax = plt.subplots(figsize = (8, 6))
    ax.errorbar(architectures, means, yerr = stds, ecolor = 'black', capsize = 10, fmt = 'o')
    print(f"Means: {means}")
    if normalize:
        print(f"Normalization factors: {distancesToSegregate}")
    ax.set_xlabel("Architectures")
    ax.set_ylim(ymin = 0)
    plt.xticks(rotation = 45)
    # ax.set_xticks(ax.get_xticks)
    # ax.set_xticklabels(ax.get_xticklabels, rotation = 45)
    if normalize:
        ax.set_ylabel("Normalized Mean Segregation Time (steps)")
        ax.set_title("Normalized Mean Segregation Times for two polymers (200 monomers each) \n of different polymer architectures")
    else:
        ax.set_ylabel("Mean Segregation Time (steps)")
        ax.set_title("Mean Segregation Times for two polymers (200 monomers each) \n of different polymer architectures")

    plt.show(block = True)
    if normalize:
        fig.savefig(f"{folder}Analysis/NormalizedSegregationTimesMeans.png")
    else:
        fig.savefig(f"{folder}Analysis/SegregationTimesMeans.png")

def AddFliersText(boxDict: dict, ax: mpl.axes.Axes) -> float:
    """Accepts the box plot dictionary returned while plotting and adds annotations to each box plot stating the number of fliers/outliers.
    Adds the annotations to the passed mpl axis."""
    for line in boxDict["fliers"]: # iterating over all box plots: each element is a Line2D object for the outliers of the box plot
        pair = line.get_data()
        numberOfFliers = len(pair[0])
        if numberOfFliers == 0:
            continue # no outliers
        x_coord = pair[0][0] # all x coords should be the same; choosing first 1
        y_coord = min(pair[1]) # choosing minimum outlier value as y coords; assuming all are upper outliers
        # Annotating the text:
        x_shift = 0.05 # a value to shift the text along the x direction
        ax.text(x_coord + x_shift, y_coord, rf"$\times{numberOfFliers}$")

def AddFliersTextAboveBoxPlots(boxDict: dict, ax: mpl.axes.Axes) -> None:
    """
    Accepts the box plot dictionary returned while plotting and adds annotations to each box plot stating the number of fliers/outliers.
    Adds the annotations to the passed mpl axis above the upper whisker of each box plot.
    """
    y_shift = 0.05
    box_props = dict(facecolor = 'white', edgecolor = 'white', alpha = 0.7, boxstyle = 'square,pad=0.01') # Textbox properties
    for i, capline in enumerate(boxDict['caps']): # iterating over all box plots
        fliers = boxDict['fliers']
        if i % 2 == 1: # Even (base 1) index corresponding to upper caps
            flierIndex = (i - 1) // 2
            numberOfFliers = len(fliers[flierIndex].get_data()[0])
            if numberOfFliers == 0:
                continue
            x_coord = fliers[flierIndex].get_data()[0][0]
            y_coord = max(capline.get_ydata())
            ax.text(x_coord, y_coord + y_shift, rf"${numberOfFliers}$", horizontalalignment = 'center', verticalalignment = 'bottom', bbox = box_props)

def GetMaxFlierCoord(boxDict: dict, areFliersAnnotated: bool) -> float:
    """Returns the y coordinate of the highest outlier/flier that is to be included in the plot"""
    y_coords = []
    for line in boxDict["fliers"]: # iterating over all box plots: each element is a Line2D object for the outliers of the box plot
        pair = line.get_data()
        numberOfFliers = len(pair[0])
        if numberOfFliers == 0:
            continue # no outliers
        if areFliersAnnotated:
            y_coord = min(pair[1]) # choosing minimum outlier (annotated) value as y coords; assuming all are upper outliers
        else:
            y_coord = max(pair[1])
        y_coords.append(y_coord)
    return max(y_coords)

def PlotBoxPlots(folder: str, architectures: list[str], ax: mpl.axis.Axis, convertToSciNotation: bool = True, showFliersInformation = True) -> tuple[dict, str]:
    """Reads and plots the distributions of segregation times for various architecutres as box plots in the passed axis
    Returns the box plot dictionary and the segregation criterion"""
    segregationTimes = [] #  a list to store segregaation times for each architecture
    segregation_criterion = "" # The criterion with which all segregations were supposedly carried out
    for architecture in architectures:
        segTimes, criterion = GetSegTimeDistribution(folder, architecture)
        if len(segregation_criterion) == 0: # not initialized
            segregation_criterion = criterion # initializing with criterion of first architecture
        else:
            if segregation_criterion != criterion: # comparing criteria of different architectures
                print(f"ERROR: All the segregation criteria for the architecures are not the same! {segregation_criterion} and {criterion} for {architecture} cannot be compared.")
                sys.exit(1)
        segTimes = np.array(segTimes)
        segTimes = segTimes / segParam.TAU_0 # converting to units in terms of TAU_0
        segregationTimes.append(segTimes)

    if convertToSciNotation:
        segregationTimes, axisLabel = pt.ConvertMultipleArraysToScientificNotation(segregationTimes, preferredOrder = segParam.PREFERRED_ORDER)
    else:
        axisLabel = ""
    # ax.violinplot(dataset=segregationTimes, showmeans=True)
    customLabels = segParam.GetAliasArchitectureList(architectures)
    xlabel = "Architecture"
    xlabel_angle = 30 # in degrees
    if segParam.USE_CUSTOM_LABELS:
        customLabels = segParam.GetCustomLabels(architectures)
        xlabel = segParam.CUSTOM_AXIS_LABEL
        xlabel_angle = 0
    boxDict = ax.boxplot(x = segregationTimes, tick_labels = customLabels)
    # Annotating outliers:
    if showFliersInformation:
        # AddFliersText(boxDict, ax)
        AddFliersTextAboveBoxPlots(boxDict, ax)
    maxFlierCoord = GetMaxFlierCoord(boxDict, showFliersInformation)

    # Setting titles and axes labels:
    simulation_label = special_simulation
    if simulation_label == "":
        simulation_label = "default"
    if segParam.SHOW_TITLE:
        ax.set_title(f"Box Plots for segregation times of different architectures\n {numberOfMonomers} monomers, {simulation_label} simulation, {segregation_criterion}")
    else:
        # ax.set_title(f"{segParam.GetAliasSimulation(special_simulation)} ({numberOfMonomers})")
        inf_label = "" # An add on label to indicate the cylinder is infinite in length
        if sysPaths.IsCylinderInfinite(special_simulation):
            inf_label = "\n" + r"$L = \infty$"
        box_props = dict(facecolor = 'white', edgecolor = 'white', alpha = 0.7)
        title_text = r"%s ($N=%i$)" % (segParam.GetAliasSimulation(special_simulation), numberOfMonomers)
        ax.text(0.6, segParam.Y_LIMS[1] * (1 - 0.2) , title_text + inf_label, fontweight = 'bold', fontsize = 26, horizontalalignment = 'left', bbox = box_props)
    ax.set_xlabel(xlabel)
    ax.tick_params(axis = 'x', labelrotation = xlabel_angle )
    ax.set_ylabel(r"$\tau_{seg}$ (%s$\tau_0$)" % axisLabel)
    # Setting y-axis limits:
    if segParam.USE_COMMON_AXIS_LIMITS:
        ax.set_ylim(bottom = segParam.Y_LIMS[0], top = segParam.Y_LIMS[1])
    else:
        up_margin = 0.1 * maxFlierCoord
        down_margin = 0.01 * maxFlierCoord
        ax.set_ylim(0 - down_margin, maxFlierCoord + up_margin)

    return boxDict, segregation_criterion

def PlotAndSaveBoxPlots(folder: str, architectures: list[str], showPlot: bool = False) -> None:
    """Plots the distributions of segregation times of all architectures passes as box plots in a single plot and saves the figure.
    The folder passed should be the location where all the architecture folders are present"""
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/big_bold.mplstyle")
    fig, ax = plt.subplots()
    boxDict, segregation_criterion = PlotBoxPlots(folder, architectures, ax, showFliersInformation = segParam.INCLUDE_FLIERS_INFORMATION)

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if segParam.USE_MARGINS:
        fig.subplots_adjust(left=segParam.PLOT_MARGINS['left'], right=segParam.PLOT_MARGINS['right'], top=segParam.PLOT_MARGINS['top'], bottom=segParam.PLOT_MARGINS['bottom'])

    # Showing and saving figure:
    if showPlot:
        plt.show(block = True)
    fig.savefig(f"{folder}{special_simulation}_{segParam.LOOP_MODIFIER}_box_segregation_times{segParam.FIG_EXT}")
    plt.close(fig)

    # Counting outliers or "fliers" in the box plots: Assuming no bottom fliers
    flierFilePath = f"{folder}{special_simulation}_fliers.csv"
    
    with open(flierFilePath, "a") as file:
        file.write("Architecture, No. of. Fliers, Flier Values\n")
        for i in range(len(boxDict["fliers"])):
            fliers = boxDict['fliers'][i].get_data()[1] # list of y-values of outliers
            file.write(f"{architectures[i]}, {len(fliers)}, {fliers}\n")

    print(f"Segregation Times Box Plots plotted for {architectures} architectures using segregation criterion {segregation_criterion}.")

def PlotBoxPlotComparison(special_simulations: list[str], architectures: list, showPlot: bool = False) -> None:
    """Plots the box plots for multiple special simulations in a single multi plot.
    Accepts a list of special simulations and a list of architectures.
    architectures can also be a 2D list of architectures with each row listing architectures for a special simulation"""

    # Handling architectures:
    arcDatabase = [[]] # A 2D List to store architectures for each special simulation
    if type(architectures[0]) == list: # CHecking if architectures is multidimensional
        arcDatabase = architectures
    else: # Duplicating architectures for each row of arcDatabase
        arcDatabase = []
        for special_simulation in special_simulations:
            arcDatabase.append(architectures)
    
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/subplots_bold.mplstyle")
    nrows = segParam.NROWS
    ncols = segParam.NCOLS
    fig, axes = plt.subplots(sharey = True, nrows = nrows, ncols = ncols)

    maxFlierCoord = 0 # To store the maximum annotated flier coord over all simulations
    hasFirstAxisBeenSet = False
    for n in range(len(special_simulations)): # Iterating over all special simulations
        folder = sysPaths.GetFolder(sysPaths.NEW_SEGREGATION, numberOfMonomers, special_simulations[n])
        i = n % nrows # row index
        j = n // nrows # column index
        if nrows == 1:
            axes_ij = axes[j]
        elif ncols == 1:
            axes_ij = axes[i]
        else:
            axes_ij = axes[i][j]
        if not hasFirstAxisBeenSet:
            axis_00 = axes_ij
            hasFirstAxisBeenSet = True
        boxDict, segregation_criterion = PlotBoxPlots(folder, arcDatabase[n], axes_ij, convertToSciNotation = False, showFliersInformation = False)
        flierCoord = GetMaxFlierCoord(boxDict, False)
        maxFlierCoord = max(maxFlierCoord, flierCoord)
        # Overwriting axis title and labels:
        axes_ij.set_title(segParam.GetAliasSimulation(special_simulations[n]))
        axes_ij.set_xlabel("")
        axes_ij.set_ylabel("")
    axis_00.set_ylim(0 - maxFlierCoord * MARGIN, maxFlierCoord + maxFlierCoord * MARGIN) # Setting the axis limits of the first subplot; teh rest should follow
    
    # Adding figure labels and title:
    fig.suptitle(f"Segregation Time Comparison, {numberOfMonomers} monomers, {segregation_criterion}")
    fig.supxlabel(segParam.CUSTOM_AXIS_LABEL)
    fig.supylabel(r"Segregation Times ($\tau_0$)")

    if showPlot:
        plt.show(block = True)
    


# script:
if __name__ == "__main__":
    SetConstants()
    directory = sysPaths.GetFolder(sysPaths.NEW_SEGREGATION, numberOfMonomers, special_simulation)
    # architectures = GetArchitectureList(directory)
    architectures = segParam.AOI_LIST
    if segParam.SORT_ARCHITECTURES:
        if segParam.USE_CUSTOM_LABELS:
            architectures.sort(key = segParam.GetCustomLabel) # sorting according to the custom labels
        else:
            architectures.sort(key = segParam.GetAliasArchitecture) # sorting according to the alias names
    # print(f"Architectures: {architectures}")
    # PlotMeans(architectures, directory, normalize = False)
    PlotAndSaveBoxPlots(directory, architectures, True)
    # PlotBoxPlotComparison(segParam.SPECIAL_SIMULATIONS, segParam.AOI_COMPARE, True)
