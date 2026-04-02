# This script is used to read a LAMMPS state file and calculate a monomer distribution
# based on just the single snapshot contained in the state file
import numpy as np
from numpy.typing import ArrayLike
import math
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from pizza import data3 # Library class to read LAMMPS state file
import sys
from pathlib import Path

# Importing additonal tools and config files
import system_file_paths as sysPaths
import Segregation_Parameters as segParam
import PlottingTools as pt
# Importing the plotMonomerDensity from another folder:
sys.path = sys.path[1:] # Excluding current directory
import plotMonomerDensity # For the AddDistributions() method to compute average
# Importing config file for regions:
sys.path.append(f"{sysPaths.POLYMER_PHYSICS}Scripts/Config_Files/")
import regions_config as reg

# Global variables:
numberOfMonomers = 200
numberOfPolymers = 2
architecture = "Arc1_10"
runIndex = -1
numberOfRuns = 50

atomsSectionHeader = "Atoms"
atomStyle = "angle"
binWidth = 0.2

plotDensity = segParam.PLOT_SINGLE_SNAPSHOT_DENSITY # A flag to indicate whether to plot the probability distribution (if True) or a histogram (if False)
diameter = 0
boxLength = 0 # Used to rescale the lengths while plotting distributions

def SetConstants():
    """Reads arguments from the command line and sets the global variables"""
    global numberOfMonomers
    global architecture
    global binWidth
    global runIndex
    global atomStyle
    global atomsSectionHeader
    global diameter
    global boxLength

    numberOfMandatoryArguments = 3
    if len(sys.argv) < numberOfMandatoryArguments + 1:
        print("Not enough arguments specified while invoking the script! Please pass the number of monomers, the name of the architecture, and the bin width in the following format:")
        print("python /<path>/plotCoMDistribution.py <noOfMonomers> <architecture> <binWidth>")
        print("Optional arguments for a run number and/or atom style can be passed to plot for just that run and a specific atomStyle.")
        print("Common atom styles include: atom, bond, angle. Default: angle")
        print(f"If 0 is passed as run number, the average over {numberOfRuns} will be performed")
        sys.exit() # terminating the script
    
    try:
        numberOfMonomers = int(sys.argv[1])
    except ValueError:
        print(f"The entered numberOfMonomers {sys.argv[1]} cannot be converted to a number! Please provide a valid number")
        sys.exit()

    architecture = sys.argv[2]

    try:
        binWidth = float(sys.argv[3])
    except ValueError:
        print(f"The entered binWidth {sys.argv[3]} cannot be converted to a number! Please provide a valid number")
        sys.exit()

    # Reading diameter and setting box length:
    diameter = pt.ReadDiameter(numberOfMonomers, architecture)
    boxLength = segParam.ASPECT_RATIO * diameter

    # Optional arguments:
    if len(sys.argv) > numberOfMandatoryArguments + 1:
        try:
            runIndex = int(sys.argv[numberOfMandatoryArguments+1])
            if len(sys.argv) > numberOfMandatoryArguments + 2:
                atomStyle = sys.argv[numberOfMandatoryArguments+2]
        except ValueError:
            if len(sys.argv) == numberOfMandatoryArguments + 2: # Just a single optional argument
                print("The optional argument is not an integer. Assuming it is the atom style.")
                atomStyle = sys.argv[numberOfMandatoryArguments+1]
            else:
                print(f"The argument {sys.argv[numberOfMandatoryArguments+1]} for the run number could not be converted to an integer! Terminating.")
                sys.exit(1)
    if len(atomStyle) > 0:
        atomsSectionHeader = f"{atomsSectionHeader} # {atomStyle}"

def GetReadFilePath(numberOfMonomers: int, architecture: str, runIndex: int) -> str:
    """Returns the file path to the mixed state created by the mixing algorithm"""
    folder = sysPaths.GetFolder(sysPaths.CREATE_INITIAL_STATES, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)
    return f"{folder}{architecture}/langevin_mixed_state_{runIndex:02}.txt"

def Read_LAMMPS_State_File(filePath: str) -> data3:
    """Reads the passed LAMMPS state file and returns a Pizza data object 
    containing the read information"""
    lmp_data = data3.data(filePath) # reading data from the state file
    # Section columns: atomID molID type x y z flag_x flag_y flag_z
    # Renaming columns for the atoms section
    lmp_data.map(1, "atomID", 2, "molID", 3, "type", 4, "x", 5, "y", 6, "z", 7, "flag_x", 8, "flag_y", 9, "flag_z")
    return lmp_data

def CalculateSingleSnapshotDistribution(data: list, binWidth: float, lower_limit: float = None, upper_limit: float = None) -> tuple[ArrayLike, ArrayLike]:
    """Calculates the probability distribution or histogram of the passed data in the range specified by lower and upper limits using the specified binwidth
    Returns the probability density and the bin_edges"""
    # Handling the case when limits are not passed:
    if sysPaths.IsCylinderInfinite():
        if lower_limit == None:
            lower_limit = min(data) - binWidth
        if upper_limit == None:
            upper_limit = max(data) + binWidth
    else:
        if lower_limit == None:
            lower_limit = -boxLength / 2
        if upper_limit == None:
            upper_limit = boxLength / 2
    # calculating bin edges:
    # binEdges = np.arange(lower_limit - binWidth, upper_limit + binWidth, binWidth) # includes an additional bin on either side
    binEdges = np.arange(lower_limit, upper_limit + binWidth, binWidth) # No extra bins at the edges; monomer density may not go to 0 at edges
    # Debugging: Printing limits and the bin edges:
    # print(f"Limits: ({lower_limit}, {upper_limit})\nBin Edges: {binEdges}")
    probDensity, binEdges = np.histogram(data, bins = binEdges, range = (lower_limit, upper_limit), density = False) # Calculates histogram; must be normalized
    if plotDensity:
        probDensity = probDensity / binWidth / (math.pi * diameter ** 2 / 4)
    
    # # Sanity check: total should be 1 / cyl_volume if density is plotted, or total count should be equal to number of monomers if histogram is plotted
    # total = np.sum(probDensity) * binWidth
    # if plotDensity:
    #     total = total * (math.pi * diameter ** 2 / 4)
    # print(f"Total integral of the computed distribution: {total}")

    return probDensity, binEdges
    
def PlotDistribution(binEdges: ArrayLike, probDensity: ArrayLike, ax: mpl.axis.Axis, label: str = "", showPlot: bool = False, putLabels: bool = True, linestyle: str = '--') -> None:
    """Plots the probability distribution or histogram using the passed numpy arrays.
    Plots on the passed axis with the given label and titles. The axis labels and title are skipped if putLabels if False."""
    binWidth = binEdges[1] - binEdges[0]
    binCentres = binEdges[:-1] + binWidth / 2
    # Rescaling lengths with box length:
    xLabelModifier = ""
    if segParam.RESCALE_LENGTHS:
        binCentres = binCentres / boxLength
        xLabelModifier = r"/ $L$"
    ax.plot(binCentres, probDensity, linestyle = linestyle, marker = ".", label = label)
    if putLabels: 
        if segParam.SHOW_TITLE:
            ax.set_title(f"Monomer distribution of 2 {architecture} polymers\n{numberOfMonomers} monomers each, 1 snapshot, Run {runIndex}")
        ax.set_xlabel(r"$z$ %s" % (xLabelModifier if len(xLabelModifier) > 0 else r"($\sigma$)"), fontweight = 'bold')
        if plotDensity:
            ax.set_ylabel(r"$n(z$%s$)/(\pi R^2 dz)$" % (xLabelModifier), fontweight = 'bold')
        else:
            ax.set_ylabel("Monomer Count", fontweight = 'bold')
    if len(label) > 0:
        leg = ax.legend()
        # Increasing legend line thickness:
        for line in leg.get_lines():
            line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())

    if showPlot:
        plt.show(block = True)

def ComputeAndPlotTotalDistribution(runIndex: int, lmp_data: data3.data, binWidth: float, lower_limit: float = None, upper_limit: float = None, showPlot: bool = False):
    """Computes the monomer probability distribution or histogram for all monomers in the system and plots it"""
    z_coords = lmp_data.get(atomsSectionHeader, 6) # z coords are in column 6

    probDensity, binEdges = CalculateSingleSnapshotDistribution(z_coords, binWidth, lower_limit, upper_limit)
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/big_bold.mplstyle")
    fig, ax = plt.subplots()
    PlotDistribution(binEdges, probDensity, ax, showPlot = showPlot)

def ReadZPositionsRegionwise(lmp_data: data3.data) -> list[list[float]]:
    """Reads the LAMMPS state file and segregates the z coords of the monomers according to the region they belong to.
    Returns a 2D list where each element list contains z coords for a region."""
    # Constructing datasets for each region
    z_coords = [] # A 2D list of coords with each row belonging to a different region
    for i in range(reg.NUMBER_OF_REGIONS):
        z_coords.append([]) # appending empty lists; z coords for atoms belonging to a particular region will be added to these lists
    totalMonomers = lmp_data.headers["atoms"]
    atomIDs = lmp_data.get(atomsSectionHeader, 1) # list of all atomIDs
    atomTypes = lmp_data.get(atomsSectionHeader, 3) # corresponding list of atom types
    allZCoords = lmp_data.get(atomsSectionHeader, 6) # corresponding list of z coords
    for n in range(totalMonomers):
        if atomTypes[n] > numberOfPolymers:
            continue # skipping monomer types which do not belong to a polymer
        regionID = reg.GetRegionID(atomIDs[n], numberOfMonomers)
        z_coords[regionID].append(allZCoords[n])
    return z_coords

def ReadAllPositionsRegionwise(lmp_data: data3.data) -> list[list[list[float]]]:
    """
    Reads the LAMMPS state file and segregates the x, y, z coords of the monomers according to the polymer they belong to.
    Returns a 3D list where each element list contains vector coords (x, y, z) for the monomers of a polymer.
    Args:
        lmp_data: Pizza data object containing the read LAMMPS state file
    """
    positions = [] # A 3D list of coords with each row belonging to a different polymer
    for i in range(reg.NUMBER_OF_REGIONS):
        positions.append([]) # appending empty lists; x, y, z coords for atoms belonging to a particular polymer will be added to these lists
    totalMonomers = lmp_data.headers["atoms"]
    # Sanity check:
    if totalMonomers != numberOfMonomers * numberOfPolymers:
        print(f"Total number of monomers in the LAMMPS state file ({totalMonomers}) does not match the expected number ({numberOfMonomers * numberOfPolymers})! Terminating.")
        sys.exit(1)
    atomIDs = lmp_data.get(atomsSectionHeader, 1) # list of all atomIDs
    atomTypes = lmp_data.get(atomsSectionHeader, 3) # corresponding list of atom types
    allXCoords = lmp_data.get(atomsSectionHeader, 4) # corresponding list of x coords
    allYCoords = lmp_data.get(atomsSectionHeader, 5) # corresponding list of y coords
    allZCoords = lmp_data.get(atomsSectionHeader, 6) # corresponding list of z coords
    for n in range(totalMonomers):
        if atomTypes[n] > numberOfPolymers:
            continue # skipping monomer types which do not belong to a polymer
        regionID = reg.GetRegionID(atomIDs[n], numberOfMonomers)
        positions[regionID].append([allXCoords[n], allYCoords[n], allZCoords[n]])
    return positions

def ComputeRegionalDistribution(region_z_coords: list[list[float]], binWidth: float, lower_limit: float = None, upper_limit: float = None) -> tuple[list[ArrayLike], list[ArrayLike]]:
    """Computes the distributions of z coords for each region. The argument is a 2D list with each element list containing z coords for a region
    Returns two lists of numpy arrays: first is the list of distribution arrays for each region, and the second is the same for bin edges of the distribution"""
    # Computing distributions:
    probDensitiesList = []
    binEdgesList = []
    for i in range(reg.NUMBER_OF_REGIONS):
        # print(f"Region {i+1}") # Debugging
        # print(f"Number of atoms: {len(z_coords[i])}") # Debugging
        probDensity, binEdges = CalculateSingleSnapshotDistribution(region_z_coords[i], binWidth, lower_limit, upper_limit)
        # Debugging: printing probability densities
        # print(f"Prob Density: {probDensity}")
        probDensitiesList.append(probDensity)
        binEdgesList.append(binEdges)
    return probDensitiesList, binEdgesList

def PlotRegionalDistributions(binEdgesList: list[ArrayLike], probDensitiesList: list[ArrayLike], runIndex: int, title: str, saveFigPath: str, showPlot: bool = False):
    """Plots the distribution for each region separately. 
    Accepts lists of distributions and the corresponding bin edges,
    run index of the run, a file path where the plot is to be saved,
    and a flag indicating whether the plot should be displayed"""
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/subplots_big_bold.mplstyle")
    fig, ax = plt.subplots(nrows = numberOfPolymers, sharex = True, sharey = True)
    for i in range(reg.NUMBER_OF_REGIONS):
        polymerIndex = i // reg.NUMBER_OF_POLYMER_REGIONS
        PlotDistribution(binEdgesList[i], probDensitiesList[i], ax[polymerIndex], reg.REGION_LABELS[i], False, False)
    # Adding titles and axis labels
    if segParam.SHOW_TITLE:
        fig.suptitle(title)
    else: # Printing a shorter title
        fig.suptitle(segParam.GetPaperFiguresTitle(numberOfMonomers, architecture))
    if segParam.RESCALE_LENGTHS:
        fig.supxlabel(r"z / $L$")
    else:
        fig.supxlabel("Coordinate along cylinder axis (LJ Units)")
    if plotDensity:
        fig.supylabel("Probability density")
    else:
        fig.supylabel("Number of Monomers")

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)
    
    # Saving figure:
    fig.savefig(saveFigPath)
    fig.subplots_adjust(bottom=0.25)
    print(f"Plotted the Single Snapshot distribution for {architecture} Run {runIndex}!")
    plt.close(fig)

def PlotRegionalDistributionsTogether(binEdgesList: list[ArrayLike], probDensitiesList: list[ArrayLike], runIndex: int, title: str, saveFigPath: str, showPlot: bool = False):
    """Plots the distribution for each region in a single axis / figure. 
    Accepts lists of distributions and the corresponding bin edges,
    run index of the run, a file path where the plot is to be saved,
    and a flag indicating whether the plot should be displayed"""
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/big_bold.mplstyle")
    fig, ax = plt.subplots(figsize = (12, 6))
    linestyles = ['--'] * reg.NUMBER_OF_REGIONS
    if reg.USE_REGIONS:
        ax.set_prop_cycle(color = [ '#87CEEB', '#00008B', '#FF6666', '#8B0000']) 
        linestyles = ['--', '--', ':', ':']
    for i in range(reg.NUMBER_OF_REGIONS):
        polymerIndex = i // reg.NUMBER_OF_POLYMER_REGIONS
        PlotDistribution(binEdgesList[i], probDensitiesList[i], ax, reg.REGION_LABELS[i], False, True, linestyle=linestyles[i])
    # Adding titles and axis labels
    if segParam.SHOW_TITLE:
        ax.set_title(title)
    else: # Printing a shorter title
        ax.text(segParam.X_LIMS[0] * (1 - 0.1), segParam.Y_LIMS[1] * (1 - 0.1) , segParam.GetPaperFiguresTitle(numberOfMonomers, architecture), fontweight = 'bold', fontsize = 26, horizontalalignment = 'left')
        # ax.text(0, segParam.Y_LIMS[1] * (1 - 0.1) , segParam.GetAliasSimulation(), fontweight = 'bold', fontsize = 26, horizontalalignment = 'center')
    # x and y labels should be set already by PlotDistribution
    
    # Setting major and minor ticks:
    # ax.xaxis.set_major_locator(MultipleLocator(0.2))
    # ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    # ax.yaxis.set_major_locator(MultipleLocator(3))
    # ax.yaxis.set_minor_locator(MultipleLocator(1))

    # Setting axis limits if necessary:
    if segParam.USE_COMMON_AXIS_LIMITS:
        ax.set_xlim(segParam.X_LIMS)
        ax.set_ylim(segParam.Y_LIMS)
    
    # Adjusting plot margins:
    if segParam.USE_MARGINS:
        fig.subplots_adjust(top = segParam.PLOT_MARGINS['top'], right = segParam.PLOT_MARGINS['right'], bottom = segParam.PLOT_MARGINS['bottom'])

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)
    
    # Saving figure:
    fig.savefig(saveFigPath)
    print(f"Plotted the Single Snapshot distribution for {architecture} Run {runIndex}!")
    plt.close(fig)

def GetFileNamePrefix() -> str:
    """Returns the file name prefix (without technial details) of the plot to be saved"""
    filePrefix = ""
    if reg.USE_REGIONS:
        filePrefix = "single_region"
    else:
        filePrefix = "single_polymer"
    distributionType = ""
    if plotDensity:
        distributionType = "_distribution"
    else:
        distributionType = "_histogram"
    return f"{filePrefix}{distributionType}"

def ComputeAndPlotRegionalDistribution(runIndex: int, lmp_data: data3.data, binWidth: float, lower_limit: float = None, upper_limit: float = None, plotTogether: bool = False, showPlot: bool = False) -> None:
    """Computes the monomer probability distribution for each region separately and plots them.
     If the plotTogether flag is False (default), the distributions are plotted on separate graphs labelled by polymer index.
     If the flag is True, all the distributions are plotted in a single figure."""
    # Reading the positions and storing them by region:
    z_coords = ReadZPositionsRegionwise(lmp_data)

    # Computing distributions:
    probDensitiesList, binEdgesList = ComputeRegionalDistribution(z_coords, binWidth, lower_limit, upper_limit)

    # Plotting distributions:
    # Setting filepath:
    fileName = GetFileNamePrefix()
    folderPath = sysPaths.GetFolder(sysPaths.CREATE_INITIAL_STATES, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)
    folderPath = "%s%s/Analysis/SingleSnapshotDistribution/" % (folderPath, architecture)
    Path(folderPath).mkdir(parents = True, exist_ok = True) # Creating the directory with parents
    filePath = "%s%s_b%02i_r%i%s" % (folderPath, fileName, int(round(binWidth * 10)), runIndex, segParam.FIG_EXT)
    # Setting title:
    title = f"Monomer Distribution of 2 {architecture} polymers\n{numberOfMonomers} monomers each, 1 snapshot, Run {runIndex}"
    if plotTogether:
        PlotRegionalDistributionsTogether(binEdgesList, probDensitiesList, runIndex, title, filePath, showPlot)
    else:
        PlotRegionalDistributions(binEdgesList, probDensitiesList, runIndex, title, filePath, showPlot)

def PlotAverageDistribution(binWidth: float, lower_limit: float = None, upper_limit: float = None, plotTogether: bool = False, showPlot: bool = False):
    """Reads, calculates, and plots the distribution averaged over all single snapshots"""
    # Preparing list of dataframes to store regionwise prob density from all runs:
    distributionDataFramesArray = [] # A 2D list with each element list corresponding to a particular region
    # Each element list contains DataFrames for all runs
    # Initializing data frames list:
    for regionIndex in range(reg.NUMBER_OF_REGIONS):
        distributionDataFramesArray.append([]) # appending empty element list for each region
    
    # The total probability distributions are stored in the following lists:
    totalProbArray = [] # A list to store the total prob distributions for each region
    totalBinsArray = [] # A list to store the bin edges for each region's total distribution

    # Filling region lists by adding prob distributions over all runs:
    for runIndex in range(1, numberOfRuns+1):
        readFilePath = GetReadFilePath(numberOfMonomers, architecture, runIndex)
        lmp_data = Read_LAMMPS_State_File(readFilePath)
        run_coords = ReadZPositionsRegionwise(lmp_data)
        # Getting the distributions for each region for a particular run:
        probList, binsList = ComputeRegionalDistribution(run_coords, binWidth, lower_limit, upper_limit)
    
        for regionIndex in range(reg.NUMBER_OF_REGIONS):
            df = pd.DataFrame(zip(binsList[regionIndex], probList[regionIndex])) # Making the dataframe object of the distribution
            # This procedure makes the length of both lists the same, so the right most edge of the bins is lost
            # We add the right most edge after the addition of distributions is done
            distributionDataFramesArray[regionIndex].append(df)
            # Debugging:
            # if runIndex == 1:
            #     print(f"Total Probability Distributions for Region {regionIndex + 1}:\n{probList[regionIndex]}")
            #     print(f"Bin Edges for Total Probability Distibution for Region {regionIndex + 1} Run {runIndex}:\n{binsList[regionIndex]}")
            #     print(f"Data Frame for Total Probability Distibution for Region {regionIndex + 1} Run {runIndex}:\n{df}")
            #     print(f"Lengths: bins = {len(binsList[regionIndex])}, probs = {len(probList[regionIndex])}")
            
    # Adding all distributions, for one region at a time:
    for regionIndex in range(reg.NUMBER_OF_REGIONS):
        sumDF = plotMonomerDensity.AddDistributions(distributionDataFramesArray[regionIndex])
        binsList = np.array(sumDF.iloc[: , 0]) # First column
        probList = np.array(sumDF.iloc[: , 1]) # Second column
        # Total probability should have been calculated
        # Adding the right most edge to the bins list:
        binWidth = binsList[1] - binsList[0]
        binsList = np.append(binsList, binWidth)
        # Normalizing with number of runs:
        probList = probList / numberOfRuns
        # Adding to the total distribution lists:
        totalBinsArray.append(binsList)
        totalProbArray.append(probList)
        # Debugging:
        # print(f"Total Probability Distributions for Region {regionIndex + 1}:\n{probList}")
        # print(f"Bin Edges for Total Probability Distibution for Region {regionIndex + 1} Run {runIndex}:\n{binsList}")
        # print(f"Lengths: bins = {len(binsList)}, probs = {len(probList)}")


    # Plotting:
    title = f"Monomer Distribution of 2 {architecture} polymers\n{numberOfMonomers} monomers each, Averaged over {numberOfRuns} snapshots"
    fileName = GetFileNamePrefix()
    folderPath = sysPaths.GetFolder(sysPaths.CREATE_INITIAL_STATES, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)
    filePath = "%s%s/Analysis/SingleSnapshotDistribution/average_%s_b%02i%s" % (folderPath, architecture, fileName, int(round(binWidth * 10)), segParam.FIG_EXT)
    if plotTogether:
        PlotRegionalDistributionsTogether(totalBinsArray, totalProbArray, runIndex, title, filePath, showPlot)
    else:
        PlotRegionalDistributions(totalBinsArray, totalProbArray, runIndex, title, filePath, showPlot)
    print(f"Plotted average monomer density of single snapshots of {architecture} polymers over {numberOfRuns} runs.")

def GetCylinderSliceLabel(sliceIndex: int, noOfLongitudinalSlices: int) -> str:
    """Returns the label for a particular longitudinal slice of the cylinder"""
    # TODO: Make this method
    sliceLabels = [r"$-\frac{L}{2} < z \leqslant -\frac{L}{4}$", r"$-\frac{L}{4} < z \leqslant 0$", r"$0 < z \leqslant \frac{L}{4}$", r"$\frac{L}{4} < z \leqslant \frac{L}{2}$"]
    return sliceLabels[sliceIndex]

def PlotCrossSectionalDistributionOnAxis(ax: mpl.axes.Axes, all_positions: list[list[list[float]]], sliceIndex: int, noOfLongitudinalSlices: int) -> None:
    """
    Plots the cross-sectional distribution of monomers in a particular longitudinal slice of the cylinder on the passed axis.
    Args:
        ax: The axis on which the cross-sectional distribution is to be plotted
        all_positions: A 3D list with each element list containing vector coords (x, y, z) for the monomers of a polymer.
        sliceIndex: The index of the longitudinal slice to be plotted (0-indexed)
        noOfLongitudinalSlices: Total number of longitudinal slices
    """
    radius = diameter / 2
    # Extracting positions in the slice:
    sliceWidth = boxLength / noOfLongitudinalSlices
    lower_limit = -boxLength / 2 + sliceIndex * sliceWidth
    upper_limit = lower_limit + sliceWidth
    markers = ['o', '^', 's', 'D', 'P', 'X'] # Different markers for different regions
    for regionIndex in range(reg.NUMBER_OF_REGIONS):
        x_coords = []
        y_coords = []
        for pos in all_positions[regionIndex]:
            if pos[2] >= lower_limit and pos[2] < upper_limit: # Checking if the z coord is in the slice
                x_coords.append(pos[0])
                y_coords.append(pos[1])
        if segParam.RESCALE_RADIUS:
            x_coords = np.array(x_coords) / radius
            y_coords = np.array(y_coords) / radius
            labelModifier = r"/ $R$"
            axisLimit = 1
        else:
            x_coords = np.array(x_coords)
            y_coords = np.array(y_coords)
            labelModifier = r" ($\sigma$)"
            axisLimit = radius
        plotLabel = ""
        if sliceIndex == 0: # first subplot
            plotLabel = reg.REGION_LABELS[regionIndex]
        ax.scatter(x_coords, y_coords, s = 100, label = plotLabel, marker = markers[regionIndex % len(markers)])
    # Plotting axis lines and cylinder boundary:
    ax.axhline(0, color = 'grey', linestyle = '--', linewidth = 1)
    ax.axvline(0, color = 'grey', linestyle = '--', linewidth = 1)
    circle = plt.Circle((0, 0), radius = axisLimit, color = 'black', fill = False, linewidth = 1)
    ax.add_patch(circle)
    ax.set_aspect('equal', 'box')
    ax.set_xlim(-axisLimit, axisLimit)
    ax.set_ylim(-axisLimit, axisLimit)
    ax.set_title(GetCylinderSliceLabel(sliceIndex, noOfLongitudinalSlices))
    if sliceIndex % segParam.NROWS == segParam.NROWS - 1: # last row
        ax.set_xlabel(r"$x$ %s" % (labelModifier))
    if sliceIndex // segParam.NROWS == 0: # first column
        ax.set_ylabel(r"$y$ %s" % (labelModifier))
    # if sliceIndex == 0: # first subplot
        # leg = ax.legend(loc = 'lower left', markerscale = 1, bbox_to_anchor = (1.2, 1), ncol = 2)
    # Setting fig legend outside the subplots later

def PlotCrossSectionalDistribution(lmp_data: data3.data, noOfLongitudinalSlices: int, runIndex: int, showPlot: bool = False) -> None:
    """
    Plots the cross-sectional distribution of monomers in different longitudinal slices of the cylinder to different subplots.
    The position distribution in the x-y plane is plotted as a 2D scatter plot, with different colours for different regions.
    """
    # Setting up the figure:
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/subplots_big_bold.mplstyle")
    # Setting bigger axis and tick labels:
    # Updating axis ticks and labels font size:
    tickLabelSize = 20
    axisLabelSize = 24
    titleSize = 24
    legendSize = 20
    mpl.rcParams.update({'xtick.labelsize': tickLabelSize, 'ytick.labelsize': tickLabelSize, 'axes.labelsize': axisLabelSize, 'axes.titlesize': titleSize, 'legend.fontsize': legendSize})

    fig, ax = plt.subplots(nrows = segParam.NROWS, ncols = segParam.NCOLS, figsize = (9, 10)) # Almost square figure
    # Reading diameter:
    diameter = pt.ReadDiameter(numberOfMonomers, architecture)
    radius = diameter / 2
    # Setting axis limits:
    # if segParam.RESCALE_LENGTHS:
    #     ax[0].set_xlim(-0.5, 0.5)
    #     ax[0].set_ylim(-0.5, 0.5)
    #     xLabelModifier = r"/ $L$"
    # else:
    #     ax[0].set_xlim(-radius, radius)
    #     ax[0].set_ylim(-radius, radius)
    #     xLabelModifier = " (LJ Units)"
    all_positions = ReadAllPositionsRegionwise(lmp_data)
    for n in range(noOfLongitudinalSlices):
        i = n % segParam.NROWS # row index
        j = n // segParam.NROWS # column index
        if segParam.NROWS == 1:
            axes_ij = ax[j]
        elif segParam.NCOLS == 1:
            axes_ij = ax[i]
        else:
            axes_ij = ax[i][j]
        PlotCrossSectionalDistributionOnAxis(axes_ij, all_positions, n, noOfLongitudinalSlices)
    # Setting legend outside the subplots:
    fig.legend(loc = 'center', markerscale = 1, ncol = 2, bbox_to_anchor = (0.5, 0.52))
    if segParam.SHOW_TITLE:
        # Adding overall title:
        fig.suptitle(f"Cross-sectional distribution of {numberOfPolymers} {architecture} polymers\n{numberOfMonomers} monomers each, Run {runIndex}")

    # Adjusting plot margins:
    if segParam.USE_MARGINS:
        fig.subplots_adjust(top = segParam.PLOT_MARGINS['top'], right = segParam.PLOT_MARGINS['right'], bottom = segParam.PLOT_MARGINS['bottom'])

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)

    # Saving figure:
    folderPath = f"{sysPaths.GetFolder(sysPaths.CREATE_INITIAL_STATES, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)}{architecture}/Analysis/SingleSnapshotDistribution/"
    filePath = f"cross_sectional_distribution_r{runIndex}{segParam.FIG_EXT}"
    Path(folderPath).mkdir(parents = True, exist_ok = True) # Creating the directory with parents
    fig.savefig(f"{folderPath}{filePath}")
    plt.close(fig)
    print(f"Plotted the cross-sectional distribution for {architecture} Run {runIndex}!")

if __name__ == "__main__":
    SetConstants()
    # Reading diameter:
    # diameter = pt.ReadDiameter(numberOfMonomers, architecture)
    # Debugging: printing diameter read from the file
    # print(f"Read diameter: {diameter}")
    if sysPaths.IsCylinderInfinite():
        limits = (None, None)
    else:
        limits = (-diameter * 5/2, diameter * 5/2) # tuple of lower limit and upper limit
    # TODO: update the limits to read from diameters database
    plotTogether = True
    if runIndex == -1:
        for runIndex in range(1, numberOfRuns+1):
            readFilePath = GetReadFilePath(numberOfMonomers, architecture, runIndex)
            lmp_data = Read_LAMMPS_State_File(readFilePath)
            ComputeAndPlotRegionalDistribution(runIndex, lmp_data, binWidth, limits[0], limits[1], plotTogether = plotTogether)
            PlotCrossSectionalDistribution(lmp_data, 4, runIndex)
    elif runIndex == 0:
        PlotAverageDistribution(binWidth, limits[0], limits[1], plotTogether = plotTogether, showPlot = False)
    else:
        readFilePath = GetReadFilePath(numberOfMonomers, architecture, runIndex)
        lmp_data = Read_LAMMPS_State_File(readFilePath)
        ComputeAndPlotRegionalDistribution(runIndex, lmp_data, binWidth, limits[0], limits[1], plotTogether = plotTogether, showPlot = True)
        # PlotCrossSectionalDistribution(lmp_data, 4, runIndex, True)

