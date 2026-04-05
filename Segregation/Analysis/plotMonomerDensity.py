import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from decimal import Decimal
import sys
import math
from numpy.typing import ArrayLike
# importing the scripts that contains all the local paths
sys.path.append(f"../../Global_Scripts/System_File_Paths/")
import system_file_paths as sysPaths 
import Segregation_Parameters as segParam
sys.path.append(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/")
import PlottingTools as pt
# importing the regions config file:
sys.path.append(f"{sysPaths.GLOBAL_SCRIPTS}Config_Files/")
import regions_config as reg

# polymer information:
numberOfPolymers = 2
indexOfPolymer = 2
architecture = ""
numberOfRuns = 50
runIndex = -1

numberOfMonomers = 200 # each 
numberOfSteps = 4 * 10 ** 7
numberOfBins = 65

boxLength = 0 # The length of the cylindrical box; used to normalize the distance along the long axis while plotting

def SetConstants():
    """Reads the argument(s) passed while invoking the script from the command line and sets it to the global variables"""
    global architecture
    global numberOfMonomers
    global runIndex
    global boxLength

    numberOfMandatoryArguments = 2
    if len(sys.argv) < numberOfMandatoryArguments + 1:
        print("Not enough arguments specified while invoking the script! Please pass the name of the architecture and the number of monomers in the following format:")
        print("python /<path>/plotMonomerDensity.py <noOfMonomers> <architecture>")
        print("Optionally, an argument for the run number can be passed after these two arguments.")
        quit() # terminating the script
    
    architecture = sys.argv[2]
    # checking if the number of monomers and regions are valid:
    try:
        numberOfMonomers = int(sys.argv[1])
    except ValueError:
        print("The first argument could not be converted to an integer!")
        sys.exit(1)

    diameter = pt.ReadDiameter(numberOfMonomers, architecture)
    boxLength = segParam.ASPECT_RATIO * diameter

    if len(sys.argv) > numberOfMandatoryArguments + 1:
        runIndex = sys.argv[numberOfMandatoryArguments+1]
        if not runIndex.isdigit():
            print(f"The entered run number {runIndex} cannot be converted to a number! Please provide a valid number")
            quit()
        runIndex = int(runIndex)

def GetFolder(special_simulation: str = sysPaths.SPECIAL_SIMULATION) -> str:
    """Returns the folder where the segregation data is stored"""
    prefix = sysPaths.GetFolder(sysPaths.SEGREGATION, numberOfMonomers, special_simulation)
    return f"{prefix}{architecture}/"

def ShiftBinEdges(binLeftEdges):
    """Accepts the list of left edges of the bins and shifts them to return a list of the bin centres"""
    binWidth = binLeftEdges[1] - binLeftEdges[0]
    binCentres = np.array(binLeftEdges) + binWidth/2
    return binCentres

def CalculateAreaUnderCurve(distribution: ArrayLike, binWidth: float, isRadial: bool = False) -> float:
    """Calculates and returns the area under the distribution curve based on the values and the bin width."""
    area = 0
    # Area = Sum(distribution[i] * binWidth)
    if not isRadial:
        area = np.sum(distribution) * binWidth
    else:
        volumePrefactor = 2 * math.pi * binWidth # prefactor for Volume of thin cylindrical shell
        sectionalAreas = [0, 0] # Area under the curve in each section
        sectionThreshold = (len(distribution) * binWidth - 0.5) / math.sqrt(2) # rough overestimate of radial threshold separating inner and outer sections
        for i in range(len(distribution)):
            averageRadius = i * binWidth + binWidth / 2 # Average of inner and outer radius for a bin
            sliverArea = distribution[i] * volumePrefactor * averageRadius
            area += sliverArea
            if i * binWidth > sectionThreshold:
                sectionalAreas[1] += sliverArea
            else:
                sectionalAreas[0] += sliverArea
        # Printing ratio of areas in sections:
        print("Debugging: Calculating area under the radial curve:")
        print(f"Inner section: {sectionalAreas[0]}\nOuter section: {sectionalAreas[1]}")
    return area

def PlotDistribution(numberOfMonomers: int, architecture: str, filePath: str, label: str, ax: mpl.axes.Axes, boxLength: float, runIndex: int = runIndex, numberOfSteps: int = numberOfSteps/segParam.TAU_0, plotDensity: bool = True) -> float:
    """Plots the monomer distribution after reading the file at filepath. The plot is made on the passed axis for a given index of polymer.
    The label argument specifies the label for the distribution being plotted.
    Returns the bin width for the distribution. The numberOfSteps is in terms of the Langevin time constant tau_0.
    The boxLength argument is used to normalize the length scales before plotting if rescaling is required.
    An optional flag plotDensity is accepted which determines if the probability density is plotted (default) 
    or the monomer counts are plotted (if False)"""

    df = pd.read_csv(filePath)

    binsLeftEdges = df.iloc[:, 0]
    binWidth = float(binsLeftEdges[1]) - float(binsLeftEdges[0])
    distribution = np.array(df.iloc[:, 1])
    
    if not plotDensity: # Converting to monomer count if the plotDensity flag is false
        distribution = distribution * binWidth * (numberOfPolymers * numberOfMonomers) # Assuming the normalization was done with total number of monomers
        # Assuming the variables numberOfPolymers and numberOfMonomers are ported from the calling script since they share the same name
    else: # Converting the distribution to scientific notation
        distribution, yAxisLabel = pt.ConvertToScientificNotation(distribution, (-1, 3))
        if len(yAxisLabel) > 0:
            yAxisLabel = f"   [{yAxisLabel}]"

    # plotting histogram:
    binCentres = ShiftBinEdges(binsLeftEdges)
    xLabelModifier = ""
    # Rescaling the bin centres by box length
    if segParam.RESCALE_LENGTHS:
        binCentres = binCentres / boxLength
        xLabelModifier = r"/ $L$"
    ax.plot(binCentres, distribution, linestyle = "--", marker = '.', label = label)
    # plt.bar(binsLeftEdges, distribution, width = binWidth, align = 'edge', alpha = 0.5, label = f'Polymer {indexOfPolymer}')
    timeLabel = rf"{Decimal(numberOfSteps):.0E}$\tau_0$"
    if segParam.SHOW_TITLE:
        ax.set_title(f"Monomer Distribution for an {architecture} polymer\nTotal: {numberOfPolymers} polymer(s), {numberOfMonomers} monomers each, {timeLabel}, Run {runIndex}")
    else: # Printing a shorter title
        ax.set_title(segParam.GetPaperFiguresTitle(numberOfMonomers, architecture))
    ax.set_xlabel(r"$z$ %s" % (xLabelModifier))
    if plotDensity:
        ax.set_ylabel(r"$P(z$%s$)$%s" % (xLabelModifier, yAxisLabel))
    else:
        ax.set_ylabel("Monomer Count")
    leg = ax.legend()
    # Increasing legend line thickness:
    for line in leg.get_lines():
        line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())
    # Debugging: Calculating area under the curve
    # print(f"{label} Area: {CalculateAreaUnderCurve(distribution, binWidth)}")

    return binWidth

def AddDistributions(distributionDataFrames: list[pd.DataFrame]) -> pd.DataFrame:
    """Accepts a list of distributions (pandas) dataframes and adds them while respecting that the 
    individual distributions might have differing extents. It extends all distributions to a common range
    before making the addition.
    Each distribution dataframe in the list must have the first column be the left bin edges
    and the next column be the distribution values.
    
    Returns:
    -A distribution dataframe with the added distributions.
    """
    # Checking validity of dataframes: must have two columns
    binWidth = -1
    for df in distributionDataFrames:
        if df.shape[1] != 2:
            print(f"ERROR: All distribution dataframes must have two columns. One of them has {df.shape[1]} column(s).")
            sys.exit(1)
        # Checking if the binwidths are the same:
        localBinWidth = df.iloc[1, 0] - df.iloc[0, 0]
        # Debugging:
        # print(f"BinWidth: {df.iloc[1, 0]} - {df.iloc[0, 0]} = {localBinWidth}")
        if binWidth < 0: # uninitialized
            binWidth = localBinWidth
        else:
            if not math.isclose(localBinWidth, binWidth):
                print(f"ERROR: The binwidths in two of the distributions are different: ({binWidth}, {localBinWidth})")
                sys.exit(1)
    
    # Assembling common range:
    binEdgeLimits = [sys.float_info.max, sys.float_info.min] # [min, max]
    for df in distributionDataFrames:
        localMax = max(df.iloc[:, 0])
        localMin = min(df.iloc[:, 0])
        if localMax > binEdgeLimits[1]:
            binEdgeLimits[1] = localMax
        if localMin < binEdgeLimits[0]:
            binEdgeLimits[0] = localMin
    commonBinEdges = np.arange(binEdgeLimits[0], binEdgeLimits[1] + binWidth, binWidth)
    # Debugging:
    # print(f"Common bin edges: {commonBinEdges}")

    # Initializing summed distribution:
    sumDist = np.zeros(len(commonBinEdges))

    # Adding distributions:
    for df in distributionDataFrames:
        for i in range(df.shape[0]): # iterating over row indices:
            commonIndex = int((df.iloc[i, 0] - binEdgeLimits[0]) // binWidth) # calculating the index of the particular bin in the common bin array
            sumDist[commonIndex] += df.iloc[i, 1]

    return pd.DataFrame({'Bin Left Edges': commonBinEdges, 'Distribution': sumDist})

def TestAddDistributions() -> None:
    """A function used to test and debug the AddDistributions() function."""
    binEdges1 = np.arange(-5, 1, 0.5)
    distribution1, binEdges1 = np.histogram(np.random.normal(loc = -2, scale = 2, size = 100), bins = binEdges1)
    df1 = pd.DataFrame({'a': binEdges1[:-1], 'b': distribution1})
    binEdges2 = np.arange(-1, 3, 0.5)
    distribution2, binEdges2 = np.histogram(np.random.normal(loc = 1, scale = 0.8, size = 100), bins = binEdges2)
    df2 = pd.DataFrame({'a': binEdges2[:-1], 'b': distribution2})

    print(f"The two bin edges:\n{binEdges1}\n{binEdges2}")
    print(f"DataFrame 1:\n{df1}")
    print(f"DataFrame 2:\n{df2}")

    sum_df = AddDistributions([df1, df2])
    # Plotting:
    fig, ax = plt.subplots(nrows = 1, ncols = 3, figsize = (16, 8), sharey = True, sharex = True)
    dists = [df1, df2, sum_df]
    labels = ["Dist 1", "Dist 2", "Sum"]
    for i in range(len(dists)):
        ax[i].plot(dists[i].iloc[:, 0], dists[i].iloc[:, 1], marker = 'o', label = labels[i])
        ax[i].legend()
    fig.supxlabel("Variable x")
    fig.supylabel("Frequency")

    plt.show(block = True)

def PlotRadialDistribution(numberOfMonomers: int, architecture: str, filePath: str, label: str, ax: mpl.axes.Axes, boxRadius: float, runIndex: int = runIndex, numberOfSteps: int = numberOfSteps/segParam.TAU_0) -> float:
    """Plots the radial monomer probability density after reading it from the file specified at filePath.
    The plot is made on the passed axis for a given index of polymer.
    The label argument specifies the label for the distribution being plotted.
    Returns the bin width for the distribution. The numberOfSteps is in terms of the Langevin time constant tau_0.
    The boxRadius argument is used to normalize the length scales before plotting if rescaling is required."""
    df = pd.read_csv(filePath)

    binsLeftEdges = df.iloc[:, 0]
    binWidth = float(binsLeftEdges[1]) - float(binsLeftEdges[0])
    distribution = np.array(df.iloc[:, 1])

    # Converting probability density to monomer number density (number of monomers per cylindrical shell volume):
    box_length = segParam.ASPECT_RATIO * boxRadius * 2
    if sysPaths.IsCylinderInfinite():
        radialNumberDensity = distribution * numberOfMonomers * numberOfPolymers # Areal density: no of monomers per unit cross sectional area
    else:
        radialNumberDensity = distribution * numberOfMonomers * numberOfPolymers / box_length # Volume density
    
    # Sanity check: Calculating volume integral of density and sum of probability:
    # total_monomers = 0
    # total_prob = 0
    # for i in range(len(radialNumberDensity)):
    #     averageRadius = binsLeftEdges[i] + binWidth / 2 # Average of inner and outer radius for a bin
    #     shellArea = 2 * math.pi * averageRadius * binWidth # Area of cross section of thin cylindrical shell
    #     shell_monomers = radialNumberDensity[i] * shellArea
    #     if not sysPaths.IsCylinderInfinite():
    #         shell_monomers = shell_monomers * box_length
    #     total_monomers += shell_monomers
    #     total_prob += distribution[i] * shellArea
    # # print(f"Debugging: binWidth = {binWidth}, boxLength = {box_length}, shellVolume (first bin) = {2 * math.pi * (binsLeftEdges[0] + binWidth / 2) * binWidth * box_length}")
    # print(f"Debugging: Total monomers calculated from radial density: {total_monomers} (Expected: {numberOfMonomers * numberOfPolymers})")
    # print(f"Debugging: Total probability from radial distribution: {total_prob} (Expected: 1.0)")

    # Converting distribution to scientific notation:
    # distribution, yAxisLabel = pt.ConvertToScientificNotation(distribution, (-1, 3))
    radialNumberDensity, yAxisLabel = pt.ConvertToScientificNotation(radialNumberDensity, (-1, 3))
    if len(yAxisLabel) > 0:
        yAxisLabel = f"   [{yAxisLabel}]"

    # plotting histogram:
    binCentres = ShiftBinEdges(binsLeftEdges)
    xLabelModifier = ""
    # Rescaling the bin centres by box length
    if segParam.RESCALE_RADIUS:
        binCentres = binCentres / boxRadius
        xLabelModifier = r"/ $R$"
    ax.plot(binCentres, radialNumberDensity, linestyle = "--", marker = '.', label = label)
    # plt.bar(binsLeftEdges, distribution, width = binWidth, align = 'edge', alpha = 0.5, label = f'Polymer {indexOfPolymer}')
    timeLabel = rf"{Decimal(numberOfSteps):.0E}$\tau_0$"
    if segParam.SHOW_TITLE:
        ax.set_title(f"Radial Monomer Distribution for an {architecture} polymer\nTotal: {numberOfPolymers} polymer(s), {numberOfMonomers} monomers each, {timeLabel}, Run {runIndex}")
    else: # Printing a shorter title
        # ax.set_title(segParam.GetPaperFiguresTitle(numberOfMonomers, architecture))
        ax.text(0.01, segParam.Y_LIMS[1] * (1 - 0.1) , segParam.GetPaperFiguresTitle(numberOfMonomers, architecture), fontweight = 'bold', fontsize = 26, horizontalalignment = 'left')
        # ax.text(0.1, 6 * (1 - 0.1) , segParam.GetAliasSimulation(), fontweight = 'bold', fontsize = 26, horizontalalignment = 'left')
    ax.set_xlabel(r"$r$ %s" % (xLabelModifier))
    if sysPaths.IsCylinderInfinite():
        volume_factor = "2 \pi r dr"
    else:
        volume_factor = "2 \pi r L dr"
    ax.set_ylabel(r"$n(r$%s$)/(%s)$%s" % (xLabelModifier, volume_factor, yAxisLabel))
    if len(label) > 0:
        leg = ax.legend()
        # Increasing legend line thickness:
        for line in leg.get_lines():
            line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())
    # Debugging: Calculating area under the curve
    print(f"{label} Area: {CalculateAreaUnderCurve(distribution, binWidth, isRadial = True)}")

    return binWidth

def PlotAndSavePolymerDistribution(showPlot: bool = False) -> None:
    """Plots the monomer distributions of all polymers and saves the figure as an image."""
    fig, ax = plt.subplots()
    folder = GetFolder()
    for indexOfPolymer in range(1, numberOfPolymers + 1):
        # getting data from CSV file:
        filePath = f"{folder}run{runIndex}/monomer_distribution_reg{indexOfPolymer}.csv"
        # filePath = f"{folder}run{runIndex}/monomer_distribution{indexOfPolymer}_0.2.csv"
        binWidth = PlotDistribution(numberOfMonomers, architecture, filePath, f"Polymer {indexOfPolymer}", ax, boxLength)

    if showPlot:
        plt.show(block = True)
    fig.savefig(f"{folder}Analysis/MonomerDistribution/MonomerDistribution_r{runIndex}{segParam.FIG_EXT}")
    plt.close(fig) # closing the figure to save memory
    print(f"Plotted Monomer Distribution for {architecture} Run {runIndex}")

def PlotAndSaveRegionDistribution(numberOfRegions: int, showPlot: bool = False) -> None:
    """Plots the monomer distributions of all the regions in the system and saves the figure as an image."""
    fig, ax = plt.subplots()
    folder = GetFolder()
    for regionIndex in range(1, numberOfRegions + 1):
        # getting data from CSV file:
        filePath = f"{folder}run{runIndex}/monomer_distribution_reg{regionIndex}.csv"
        binWidth = PlotDistribution(numberOfMonomers, architecture, filePath, reg.REGION_LABELS[regionIndex-1], ax, boxLength)

    if showPlot:
        plt.show(block = True)
    fig.savefig(f"{folder}Analysis/MonomerDistribution/MonomerDistribution_r{runIndex}{segParam.FIG_EXT}")
    plt.close(fig) # closing the figure to save memory
    print(f"Plotted Monomer Distribution for {architecture} Run {runIndex}")

def PlotAndSaveRadialDistribution(showPlot: bool = False) -> None:
    """Plots the radial monomer distribution and saves the figure as an image"""
    folder = GetFolder()
    filePath = f"{folder}run{runIndex}/radial_monomer_distribution.csv"
    radius = boxLength / segParam.ASPECT_RATIO / 2
    fig, ax = plt.subplots()
    binWidth = PlotRadialDistribution(numberOfMonomers, architecture, filePath, "", ax, radius, runIndex)
    
    if showPlot:
        plt.show(block = True)
    fig.savefig(f"{folder}Analysis/RadialDistribution/RadialDistribution_r{runIndex}{segParam.FIG_EXT}")
    plt.close(fig) # closing the figure to save memory
    print(f"Plotted Radial Monomer Distribution for {architecture} Run {runIndex}")

def PlotRadialInitializationComparison(special_simulations: list[str], numberOfMonomers: int, architecture: str, boxLength: float, runIndex: int = runIndex, baseFolder: str = sysPaths.SEGREGATION, showPlot: bool = False) -> None:
    """
    Plots a comparison of the radial monomer distributions for different initialization conditions of the same architecture system.
    Args:
        - special_simulations: A list of strings specifying the special simulation types (e.g., ["fene_recenter", "fene_fbsr"])
        - numberOfMonomers: Number of monomers in each polymer
        - architecture: The architecture of the polymer (string)
        - boxLength: The (finite) length of the cylindrical box
        - runIndex: The index of the run to be plotted. Default is the global variable runIndex.
        - baseFolder: The base folder where the run data is stored. Default is sysPaths.SEGREGATION
        - showPlot: A boolean flag to indicate if the plot should be displayed interactively. Default is False.
    """
    fig, ax = plt.subplots()
    boxRadius = boxLength / segParam.ASPECT_RATIO / 2
    subfolder = ""
    if baseFolder == sysPaths.CREATE_INITIAL_STATES:
        subfolder = "monomer_distribution/"
    for special_simulation in special_simulations:
        folder =f"{sysPaths.GetFolder(baseFolder, numberOfMonomers, special_simulation)}{architecture}/"
        filePath = f"{folder}run{runIndex}/{subfolder}radial_monomer_distribution.csv"
        df = pd.read_csv(filePath)

        binsLeftEdges = df.iloc[:, 0]
        binWidth = float(binsLeftEdges[1]) - float(binsLeftEdges[0])
        distribution = np.array(df.iloc[:, 1])

        # Converting probability density to monomer number density (number of monomers per cylindrical shell volume):
        if sysPaths.IsCylinderInfinite():
            radialNumberDensity = distribution * numberOfMonomers * numberOfPolymers # Areal density: no of monomers per unit cross sectional area
        else:
            radialNumberDensity = distribution * numberOfMonomers * numberOfPolymers / boxLength # Volume density
        
        # Converting distribution to scientific notation:
        radialNumberDensity, yAxisLabel = pt.ConvertToScientificNotation(radialNumberDensity, (-1, 3))
        if len(yAxisLabel) > 0:
            yAxisLabel = f"   [{yAxisLabel}]"

        # plotting histogram:
        binCentres = ShiftBinEdges(binsLeftEdges)
        xLabelModifier = ""
        # Rescaling the bin centres by box length
        if segParam.RESCALE_RADIUS:
            binCentres = binCentres / boxRadius
            xLabelModifier = "/ R"
        ax.plot(binCentres, radialNumberDensity, linestyle = "--", marker = '.', label = segParam.GetAliasSimulation(special_simulation))
        # plt.bar(binsLeftEdges, distribution, width = binWidth, align = 'edge', alpha = 0.5, label = f'Polymer {indexOfPolymer}')
        timeLabel = rf"{Decimal(numberOfSteps):.0E}$\tau_0$"
        if segParam.SHOW_TITLE:
            ax.set_title(f"Radial Monomer Distribution for an {architecture} polymer\nTotal: {numberOfPolymers} polymer(s), {numberOfMonomers} monomers each, {timeLabel}, Run {runIndex}")
        ax.set_xlabel(r"$r$ %s" % (xLabelModifier))
        if sysPaths.IsCylinderInfinite():
            volume_factor = "2 \pi r dr"
        else:
            volume_factor = "2 \pi r L dr"
        # ax.set_ylabel(r"$\left\langle n(r%s) \right\rangle \: / \: 2N$%s" % (xLabelModifier, yAxisLabel))
        ax.set_ylabel(r"$n(r%s) / (%s)$ %s" % (xLabelModifier, volume_factor, yAxisLabel))
        # Debugging: Calculating area under the curve
        # print(f"{label} Area: {CalculateAreaUnderCurve(distribution, binWidth, isRadial = True)}")
    leg = ax.legend()
    # Setting axis limits if necessary:
    if segParam.USE_COMMON_AXIS_LIMITS:
        ax.set_ylim(top = 6)

    # Increasing legend line thickness:
    for line in leg.get_lines():
        line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())

    # Adjusting plot margins:
    if segParam.USE_MARGINS:
        fig.subplots_adjust(top = segParam.PLOT_MARGINS['top'], right = segParam.PLOT_MARGINS['right'], bottom = segParam.PLOT_MARGINS['bottom'])

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)
    
    folder = f"{sysPaths.GetFolder(baseFolder, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)}{architecture}/"
    fig.savefig(f"{folder}Analysis/RadialDistribution/RadialDistribution_InitializationComparison_r{runIndex}{segParam.FIG_EXT}")
    plt.close(fig) # closing the figure to save memory
    print(f"Plotted Radial Monomer Distribution Comparison for {architecture} Run {runIndex}")

def PlotAndSaveRegionwiseRadialDistribution(showPlot: bool = False) -> None:
    """Plots the radial monomer distributions of all the regions in the system and saves the figure as an image."""
    fig, ax = plt.subplots()
    folder = GetFolder()
    for regionIndex in range(1, reg.NUMBER_OF_REGIONS + 1):
        # getting data from CSV file:
        filePath = f"{folder}run{runIndex}/radial_monomer_distribution_reg{regionIndex}.csv"
        radius = boxLength / segParam.ASPECT_RATIO / 2
        binWidth = PlotRadialDistribution(numberOfMonomers, architecture, filePath, reg.REGION_LABELS[regionIndex-1], ax, radius, runIndex)

    if showPlot:
        plt.show(block = True)
    fig.savefig(f"{folder}Analysis/RadialDistribution/regionwiseDistribution_r{runIndex}.png")
    plt.close(fig) # closing the figure to save memory
    print(f"Plotted regionwise Radial Monomer Distribution for {architecture} Run {runIndex}")

#script:
if __name__ == "__main__":
    SetConstants()
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/bold.mplstyle")
    if runIndex == -1:
        for runIndex in range(1, numberOfRuns + 1):
            if reg.USE_REGIONS:
                PlotAndSaveRegionDistribution(reg.NUMBER_OF_REGIONS)
            else:
                PlotAndSavePolymerDistribution()
            if segParam.PLOT_COMPARISON:
                PlotRadialInitializationComparison(segParam.SPECIAL_SIMULATIONS, numberOfMonomers, architecture, boxLength, runIndex, True)
            else:
                PlotAndSaveRadialDistribution()
        print(f"Successfully plotted the Monomer Distributions for the {numberOfRuns} runs of {architecture}")
    else:
        if reg.USE_REGIONS:
            PlotAndSaveRegionDistribution(reg.NUMBER_OF_REGIONS, True)
        else:
            PlotAndSavePolymerDistribution(showPlot = True)
        if segParam.PLOT_COMPARISON:
            PlotRadialInitializationComparison(segParam.SPECIAL_SIMULATIONS, numberOfMonomers, architecture, boxLength, runIndex, True)
        else:
            PlotAndSaveRadialDistribution(True)
        PlotAndSaveRegionwiseRadialDistribution(True)
        print(f"Succesfully plotted the Monomer Distributions for the Run {runIndex} of {architecture}")
