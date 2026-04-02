import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
from pathlib import Path
# importing the scripts that contains all the local paths
# The path to the directory containing this module is set as the environment variable PYTHONPATH
import system_file_paths as sysPaths 
import Segregation_Parameters as segParam
# Importing the plotMonomerDensity from New_Segregation/Analysis directory:
sys.path = sys.path[1:] # Excluding current directory
import plotMonomerDensity
# importing PlottingTools:
# The path to the directory containing this module is set as the environment variable PYTHONPATH
import PlottingTools as pt

# Global variables:
numberOfMonomers = 200
architecture = "Arc0"
runIndex = -1
binWidth = 0.1
plotPairCorrelationSeparately = False

forSegregation = False #  A flag to indicate whether the pair correlation should be plotted for a segregation run

acceptedModes = ["monomer", "loop"] # A list of the accepted modes that the script can operate in
chosenMode = ""

numberOfRuns = 50

def SetConstants() -> None:
    """Reads the arguments passed to the command line and sets the values of the global variables"""
    global numberOfMonomers
    global architecture
    global runIndex
    global chosenMode
    global binWidth
    global plotPairCorrelationSeparately
    global forSegregation

    numberOfMandatoryArguments = 6
    if len(sys.argv) < numberOfMandatoryArguments + 1:
        print("Not enough arguments specified while invoking the script! Please pass the number of monomers, the name of the architecture, the bin width, the mode of operation, the flag to plot the pair correlation separately, and the flag for segregation in the following format:")
        print(f"The operation modes should be chosen from the following options: {acceptedModes}")
        print("python /<path>/plotPairCorrelation.py <noOfMonomers> <architecture> <binWidth> <mode> <plot-flag> <segregation-flag>")
        print("If the plot-flag is true, the following is plotted:\n-pair correlation between monomers of the same polymer,\n-pair correlation between monomers from different polymers, and\n-total pair correlation function.")
        print("If the segregation-flag is true, then the pair correlation is plotted for a segregation run. Otherwise, it is plotted for the initial state run.")
        print("An optional argument for a run number can be passed to plot for just that run.")
        sys.exit(1)

    try:
        numberOfMonomers = int(sys.argv[1])
        binWidth = float(sys.argv[3])
    except ValueError:
        print(f"One (or both) of the two arguments: {sys.argv[1]} and {sys.argv[3]} could not be converted to numbers! Terminating.")
        sys.exit(1)
    
    architecture = sys.argv[2]

    chosenMode = sys.argv[4]
    if chosenMode not in acceptedModes:
        print(f"ERROR: The mode '{chosenMode}' was not recognised. Please pass a mode from one of the following: {acceptedModes}")
        sys.exit(1)

    if sys.argv[5] == "True" or sys.argv[5] == "1":
        plotPairCorrelationSeparately = True
    elif sys.argv[5] == "False" or sys.argv[5] == "0":
        plotPairCorrelationSeparately = False
    else:
        print(f"The passed argument {sys.argv[5]} for the flag could not be interpreted as a boolean value!")
        sys.exit(1)

    if sys.argv[6] == "True" or sys.argv[6] == "1":
        forSegregation = True
    elif sys.argv[6] == "False" or sys.argv[6] == "0":
        forSegregation = False
    else:
        print(f"The passed argument {sys.argv[6]} for the flag could not be interpreted as a boolean value!")
        sys.exit(1)

    # Optional argument
    if len(sys.argv) > numberOfMandatoryArguments + 1 :
        try:
            runIndex = int(sys.argv[numberOfMandatoryArguments+1])
        except ValueError:
            print(f"The optional argument for run index {sys.argv[numberOfMandatoryArguments+1]} could not be converted to an integer! Terminating.")
            sys.exit(1)

def GetFolder(baseFolder: str, numberOfMonomers: int, architecture: str) -> str:
    return f"{sysPaths.GetFolder(baseFolder, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)}{architecture}/"

def GetPairCorrelationFilePath(folder: str, chosenMode: str, runIndex: int, binWidth: float) -> str:
    """
    Returns the file path to the file where the pair correlation function is written.
    Args:
        -folder: The path to the directory containing the data for all simulaitons with a given number of monomers (For eg: path to the 'b200/' directory)
        -chosenMode: One of the string modes that is accepted by the script from the command line
        -runIndex: The index of the run folder for which the pair corrrelation is to be plotted
        -binWidth: The bin width with which the pair correlation was calculated.
    """
    if chosenMode == "monomer":
        readFilePath = f"{folder}run{runIndex}/pair_correlation/pair_correlation_b{int(round(binWidth*10)):02}.csv"
    elif chosenMode == "loop":
        readFilePath = f"{folder}run{runIndex}/pair_correlation/loops_pair_correlation_b{int(round(binWidth*10)):02}.csv"
    else:
        print(f"The chosen mode '{chosenMode}' could not be recognised. Please pass a mode from the following: {acceptedModes}")
        sys.exit(1)
    return readFilePath

def ReadAndPlotPairCorrelation(folder: str, runIndex: int, binWidth: float, plotSeparately: bool = False, showPlot: bool = False) -> None:
    """Reads the calculated pair correlation function from the given folder and plots it.
    Accepts a runIndex, binwidth as positional arguments.
    Accepts a flag as a keyword argument to plot the pair correlation separately: 
    -pair correlation between monomers of the same polymer, and
    -pair correlation between monomers from different polymers.
    Accepts a flag to show the plotted figure as a keyword argument."""
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/big_bold.mplstyle")
    fileNames = ["pair_correlation"] # only one file present by default
    labels = ["Total"]
    if plotSeparately:
        fileNames.append("same_pair_correlation")
        fileNames.append("diff_pair_correlation")
        labels = ["Total", "Same", "Different"]

    fig, ax = plt.subplots()

    # Calculating the sum of the two separate pair correlation functions; for debugging
    # sum = []
    # isSumInitialized = False # A flag to indicate whether the sum array has been initialized
    for i in range(len(fileNames)):
        df = pd.read_csv("%srun%i/pair_correlation/%s_b%02i.csv" % (folder, runIndex, fileNames[i], binWidth * 10))
        leftBinEdges = df.iloc[: , 0]
        pairCorrelation = df.iloc[: , 1]
        pairCorrelation = np.array(pairCorrelation) * 2
        # if i > 0: # Adding to sum only if the pair correlation is one of the same or diff
        #     if not isSumInitialized:
        #         sum = np.array(pairCorrelation) # initializing with the first set of pair correlation
        #         isSumInitialized = True
        #     else:
        #         sum += np.array(pairCorrelation)
        calcbinWidth = leftBinEdges[1] - leftBinEdges[0] # calculated binWidth
        if binWidth != calcbinWidth:
            print("The passed bin width does not match the bin width in the data in the %s file!" % (fileNames[i]))
        binCentres = plotMonomerDensity.ShiftBinEdges(leftBinEdges)
        # ax.bar(leftBinEdges, pairCorrelation, binWidth, align = 'edge')
        ax.plot(binCentres, pairCorrelation, linestyle = '--', marker = '.', label = labels[i])
    # Plotting sum pair correlation: debugging
    # ax.plot(binCentres, sum , linestyle = '--', marker = '.', label = "Sum")
    if segParam.SHOW_TITLE:
        ax.set_title(f"Pair Correlation function for two {architecture} polymers\n{numberOfMonomers} monomers each, Bin Width = {binWidth}, Run {runIndex}")
    else: # Printing a shorter title
        ax.set_title(segParam.GetPaperFiguresTitle(numberOfMonomers, architecture))
    ax.set_ylabel("g(r)")
    ax.set_xlabel("Distance r")
    leg = ax.legend()
    for line in leg.get_lines():
        line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())
    fig.subplots_adjust(bottom=0.15)
    
    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)

    # Saving the figure:
    saveFolder = "%sAnalysis/PairCorrelation/" % (folder)
    Path(saveFolder).mkdir(parents = True, exist_ok = True) # Making the directory
    fig.savefig("%spair_correlation_b%02i_r%i%s" % (saveFolder, int(round(binWidth*10)), runIndex, segParam.FIG_EXT)) # Multiply 10 due to a convention
    plt.close(fig)
    print(f"Plotted pair correlation function for {architecture} Run {runIndex}.")

def ReadAndPlotDistanceDistribution(folder: str, runIndex: int, binWidth: float, showPlot: bool = False) -> None:
    """Reads the calculated distance distribution from the given folder and plots it"""

    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/bold.mplstyle")

    df = pd.read_csv("%srun%i/pair_correlation/distance_distribution_b%02i.csv" % (folder, runIndex, binWidth * 10))
    leftBinEdges = df.iloc[: , 0]
    distanceDistribution = df.iloc[: , 1]
    calcbinWidth = leftBinEdges[1] - leftBinEdges[0] # calculated binWidth
    if binWidth != calcbinWidth:
        print("The passed bin width does not match the bin width in the data!")
    fig, ax = plt.subplots()
    binCentres = plotMonomerDensity.ShiftBinEdges(leftBinEdges)
    # print(binCentres)
    # ax.bar(leftBinEdges, pairCorrelation, binWidth, align = 'edge')
    ax.plot(binCentres, distanceDistribution, linestyle = '--', marker = '.')

    if segParam.SHOW_TITLE:
        ax.set_title(f"Pairwise Distance distribution for two {architecture} polymers\n{numberOfMonomers} monomers each, Bin Width = {binWidth}, Run {runIndex}")
    ax.set_ylabel("Probability distribution")
    ax.set_xlabel("Distance r")

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)

    # Saving the figure:
    saveFolder = "%sAnalysis/PairCorrelation/" % (folder)
    Path(saveFolder).mkdir(parents = True, exist_ok = True) # Making the directory
    fig.savefig("%sdist_distribution_b%02i_r%i%s" % (saveFolder, int(round(binWidth*10)), runIndex, segParam.FIG_EXT))
    plt.close(fig)
    print(f"Plotted distance distribution for {architecture} Run {runIndex}.")

def ReadAndPlotLoopsPairCorrelation(folder: str, runIndex: int, binWidth: float, showPlot: bool = False) -> None:
    """
    Reads the pair corretion function written for the small subloops as particles and plots it.
    """
    # Setting file path:
    readFilePath = GetPairCorrelationFilePath(folder, chosenMode, runIndex, binWidth)
    # Reading file:
    df = pd.read_csv(readFilePath)
    binLeftEdges = df.iloc[: , 0]
    pairCorrelation = df.iloc[: , 1]
    binCentres = plotMonomerDensity.ShiftBinEdges(binLeftEdges)

    # Plotting:
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/bold.mplstyle")
    fig, ax = plt.subplots()
    ax.plot(binCentres, pairCorrelation, linestyle = '--', marker = '.')
    ax.set_title(f"Pair Correlation function for subloops of two {architecture} polymers\n{numberOfMonomers} monomers each, Bin Width = {binWidth}, Run {runIndex}")
    ax.set_ylabel("Pair Correlation Function g(r)")
    ax.set_xlabel("Distance r")

    if showPlot:
        plt.show(block = True)
    # Saving figure:
    saveFolder = "%sAnalysis/PairCorrelation/" % (folder)
    Path(saveFolder).mkdir(parents = True, exist_ok = True) # Making the directory
    fig.savefig("%sloops_pair_correlation_b%02i_r%i%s" % (saveFolder, int(round(binWidth*10)), runIndex, segParam.FIG_EXT)) # Multiply 10 due to a convention
    plt.close(fig)
    print(f"Plotted loops pair correlation function for {architecture} Run {runIndex}.")

def PlotComparison(runIndex: int, showPlot: bool = False) -> None:
    """
    Plots the pair correlation functions for multiple architectures given in AOI_COMPARE defined Segregation_Parameters.py.
    Plots the function in a single figure.
    Args:
        -runIndex: The run index for which the pair correlation functions are read for all architectures
    """
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/bold.mplstyle")
    fig, ax = plt.subplots()
    if forSegregation:
        baseFolder = sysPaths.NEW_SEGREGATION
    else:
        baseFolder = sysPaths.CREATE_INITIAL_STATES

    for arc in segParam.AOI_COMPARE:
        folder = GetFolder(baseFolder, numberOfMonomers, arc)
        readFilePath = GetPairCorrelationFilePath(folder, chosenMode, runIndex, binWidth)
        df = pd.read_csv(readFilePath)
        binLeftEdges = df.iloc[: , 0]
        pairCorrelation = df.iloc[: , 1]
        binCentres = plotMonomerDensity.ShiftBinEdges(binLeftEdges)

        # Plotting:
        ax.plot(binCentres, pairCorrelation, linestyle = '--', marker = '.', label = f"{arc}")
    # Setting titles and labels:
    loopStr = ""
    if chosenMode == "loop":
        loopStr = "for subloops "
    ax.set_title(f"Pair Correlation function {loopStr}of two polymers\n{numberOfMonomers} monomers each, Bin Width = {binWidth}, Run {runIndex}")
    ax.set_ylabel("Pair Correlation Function g(r)")
    ax.set_xlabel("Distance r")
    plt.legend()

    if showPlot:
        plt.show(block = True)
    
    # Saving figure:
    folder = GetFolder(baseFolder, numberOfMonomers, architecture)
    saveFolder = "%sAnalysis/PairCorrelation/" % (folder)
    Path(saveFolder).mkdir(parents = True, exist_ok = True) # Making the directory
    fig.savefig("%scomparison_pair_correlation_b%02i_r%i%s" % (saveFolder, int(round(binWidth*10)), runIndex, segParam.FIG_EXT)) # Multiply 10 due to a convention
    plt.close(fig)
    print(f"Plotted comparison of pair correlation functions for {segParam.AOI_COMPARE} Run {runIndex}.")

    


if __name__ == "__main__":
    SetConstants()
    # Testing SetConstants:
    # print(f"Number of Monomers: {numberOfMonomers}")
    # print(f"Architecture: {architecture}")
    # print(f"Binwidth: {binWidth}")
    # print(f"Chosen Mode: {chosenMode}")
    # print(f"Plotting flag: {plotPairCorrelationSeparately}")
    # print(f"Segregation flag: {forSegregation}")
    # if runIndex != -1:
    #     print(f"Run Index: {runIndex}")
    folder = GetFolder(sysPaths.CREATE_INITIAL_STATES, numberOfMonomers, architecture)
    if forSegregation:
        folder = GetFolder(sysPaths.NEW_SEGREGATION, numberOfMonomers, architecture)
    # print(f"Folder: {folder}")
        
    if runIndex == -1:
        for i in range(1, numberOfRuns+1):
            if segParam.PLOT_COMPARISON:
                PlotComparison(i)
            else:
                if chosenMode == "monomer":
                    ReadAndPlotPairCorrelation(folder, i, binWidth, plotSeparately = plotPairCorrelationSeparately)
                else: # chosenMode == "loop"
                    ReadAndPlotLoopsPairCorrelation(folder, i, binWidth)
        print(f"Plotted pair correlation functions for {numberOfRuns} runs of {architecture}.")
    else: # Optional run index passed:
        if segParam.PLOT_COMPARISON:
            PlotComparison(runIndex, showPlot = True)
        else:
            if chosenMode == "monomer":
                ReadAndPlotPairCorrelation(folder, runIndex, binWidth, plotSeparately = plotPairCorrelationSeparately, showPlot = True)
            else: # chosenMode == "loop"
                ReadAndPlotLoopsPairCorrelation(folder, runIndex, binWidth, showPlot = True)