import pandas as pd
import numpy as np
from numpy.typing import ArrayLike
import matplotlib as mpl
import matplotlib.pyplot as plt
from decimal import Decimal
import sys
import math
from pathlib import Path
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator) # For setting major and minor ticks
# importing the scripts that contains all the local paths
# The path to the directory containing this module is set as the environment variable PYTHONPATH
sys.path.append(f"../../Global_Scripts/System_File_Paths/") # Adding the path to the system file paths module to the system path
import system_file_paths as sysPaths 
sys.path.append(f"{sysPaths.SEGREGATION}Analysis/")
import Segregation_Parameters as segParam
sys.path.append(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/")
import PlottingTools as pt
# Importing the other single snapshot distribution script
import plotSingleSnapshotDistribution as single
# Importing the plotMonomerDensity from another folder:
sys.path = sys.path[1:] # Excluding current directory
import plotMonomerDensity
sys.path.append(f"{sysPaths.GLOBAL_SCRIPTS}Config_Files/")
import regions_config as reg

# polymer information:
numberOfPolymers = 2
indexOfPolymer = 2
architecture = ""
numberOfRuns = 50
runIndex = -1

numberOfSteps = 4 * 10 ** 7
TAU_0 = 200 # timestep = 0.05 and damp = 1

initializationProcedure = ""

numberOfMonomers = 200 # each 
diameter = 0
boxLength = 0

def SetArguments():
    """Reads the argument(s) passed while invoking the script from the command line and sets it to the global variables"""
    global architecture
    global numberOfMonomers
    global runIndex
    global diameter
    global boxLength
    global initializationProcedure
    
    numberOfMandatoryArguments = 3
    if len(sys.argv) < numberOfMandatoryArguments + 1:
        print("Not enough arguments specified while invoking the script! Please pass the number of monomers, the name of the architecture, and the name of the initialization procedure in the following format:")
        print("python /<path>/plotMonomerDensity.py <noOfMonomers> <architecture> <initializationProcedure>")
        print("An optional argument for a run number can be passed to plot for just that run.")
        print("If 0 is passed as run index, the average monomer distribution over all runs will be plotted (only when segParam.USE_SINGLE_SNAPSHOT = True).")
        sys.exit(1) # terminating the script
    
    architecture = sys.argv[2]
    initializationProcedure = sys.argv[3]
    numberOfMonomers = sys.argv[1]

    # checking if the number of monomers is valid:
    if not numberOfMonomers.isdigit():
        print(f"The entered numberOfMonomers {numberOfMonomers} cannot be converted to a number! Please provide a valid number")
        sys.exit(1)

    numberOfMonomers = int(numberOfMonomers)
    
    # Reading Box Length:
    diameter = pt.ReadDiameter(numberOfMonomers, architecture, initializationProcedure)
    boxLength = segParam.ASPECT_RATIO * diameter

    # Optional argument:
    if len(sys.argv) > numberOfMandatoryArguments + 1:
        try:
            runIndex = int(sys.argv[numberOfMandatoryArguments+1])
        except ValueError:
            print("The argument for the run number could not be converted to an integer! Terminating.")
            sys.exit(1)


def ShiftBinEdges(binLeftEdges, binWidth: float):
    """Takes in the collection of bin edges and shifts to make them the sin centres"""
    binEdges = np.array(binLeftEdges)
    binCentres = binEdges + 0.5 * binWidth
    return binCentres

def GetFolder(numberOfMonomers: int, architecture: str, special_simulation: str) -> str:
    """Returns the path to the folder where the distribution files are stored"""
    folderPrefix = sysPaths.GetFolder(sysPaths.CREATE_INITIAL_STATES, numberOfMonomers, special_simulation)
    folder = f"{folderPrefix}{architecture}/"
    return folder

def GetDistributionFilePath(folder: str, runIndex: int, fileName: str) -> str:
    """Returns the the path to the monomer distribution file with the specifed run and atom type"""
    # return f"{folder}run{runIndex}/monomer_distribution{fileIndex}.csv"
    return f"{folder}run{runIndex}/monomer_distribution/{fileName}"

def GetDistributionData(filePath: str) -> tuple[ArrayLike, ArrayLike]:
    """Reads the distribution file and returns a tuple of the bin centres and corresponding distribution values"""
    df = pd.read_csv(filePath)
    binsLeftEdges = df.iloc[:, 0]
    binWidth = float(binsLeftEdges[1]) - float(binsLeftEdges[0])
    # shifting the binsLeftEdges so that they become the midpoint of the bins
    binCentres = ShiftBinEdges(binsLeftEdges, binWidth)
    distribution = df.iloc[:, 1]
    return binCentres, distribution

def PlotSingleDistribution(ax, binCentres: ArrayLike, distribution: ArrayLike, title: str) -> None:
    """Plots the passed data on the passed axis"""
    ax.plot(binCentres, distribution, linestyle = '--', marker = '.')
    if len(title) > 0:
        ax.set_title(title)
    # ax.set_xlabel("z coordinate along cylinder axis (simulation units)")
    # ax.set_ylabel("Probability Density")

    # debugging the skewness of the data set:
    # mean = np.mean(distribution)
    # ax.axhline(y = mean, linestyle = ':', color = 'r')
    
    ax.tick_params(axis = 'both', labelleft = True, labelbottom = True) # ensuring that the usual x and y tick labels are visible

def PlotDistribution(showPlot: bool = False) -> None:
    """Plots the monomer distribution for an entire system"""
    fig, ax = plt.subplots(figsize = (10, 6))
    folder = GetFolder(numberOfMonomers, architecture, special_simulation=initializationProcedure)
    filePath = GetDistributionFilePath(folder, runIndex, f"radial_monomer_distribution.csv")
    binCentres, distribution = GetDistributionData(filePath)
    PlotSingleDistribution(ax, binCentres, distribution, "")
    binWidth = binCentres[1] - binCentres[0]
    ax.set_title("Monomer Distribution for %i %s polymer(s)\n with %i each, Run %i, bin width = %.1lf" % (numberOfPolymers, architecture, numberOfMonomers, runIndex, binWidth))
    # fig.tight_layout()

    if showPlot:
        plt.show(block = True)
    # fig.savefig(f"{folder}run{runIndex}/monomer_dist.png")

def ConvertToMonomerNumberDensity(distribution: ArrayLike) -> ArrayLike:
    """Converts the passed probability distribution for a polymer/region/section to monomer number density"""
    numberDensity = np.array(distribution) * (numberOfMonomers * numberOfPolymers / (math.pi * diameter ** 2 / 4))
    return numberDensity

def PlotDistributionForAllSections(numberOfSections: int, showPlot: bool = False) -> None:
    """Plots the distributions for all sections in a cylinder"""
    nrows = 2
    ncols = numberOfSections // nrows
    fig, axes = plt.subplots(nrows = nrows, ncols = ncols, sharey = True, sharex = True)
    sectionID = 1
    folder = GetFolder(numberOfMonomers, architecture, special_simulation=initializationProcedure)
    # debugging:
    # totalProbability = 0 # verifying that that the area under the prob distribution curves amounts to 1
    # binWidth = 0.4
    distribution_list = []
    number_density_list = []
    binCentres_list = []
    vol_integral = 0
    for i in range(nrows):
        for j in range(ncols):
            filePath = GetDistributionFilePath(folder, runIndex, f"monomer_distribution_sect{sectionID}.csv")
            binCentres, distribution = GetDistributionData(filePath)

            # Converting to monomer number density:
            number_density = ConvertToMonomerNumberDensity(np.array(distribution))
            # Debugging:
            binCentres_list.append(binCentres)
            distribution_list.append(distribution) # Adding all distributions in a list to convert to scientific notation
            number_density_list.append(number_density)
            sectionID += 1

    # # Sanity check: Ensuring that the volume integral of all number densities add up to the appropriate number of monomers
    # binWidth = binCentres[1] - binCentres[0]
    # total_monomers = 0
    # for indexOfPSection in range(len(number_density_list)):
    #     section_monomers = np.sum(number_density_list[indexOfPSection]) * binWidth * (math.pi * diameter ** 2 / 4)
    #     print(f"Total monomers in Section {indexOfPSection + 1}: {section_monomers}")
    #     total_monomers += section_monomers
    # print(f"Total monomers in all sections: {total_monomers} (Should be close to {numberOfMonomers * numberOfPolymers})")

    # distribution_list, axisLabel = pt.ConvertMultipleArraysToScientificNotation(distribution_list, (-1, 3))
    number_density_list, axisLabel = pt.ConvertMultipleArraysToScientificNotation(number_density_list, (-1, 3))
    sectionID = 1 # Counter to keep track of the index/ sectionID of the distributions
    binWidth = binCentres[1] - binCentres[0]
    # Debugging: Printing binwidth
    # print("Sectional Density Binwidth: %.12lf, binWidth label: %02i" % (binWidth, binWidth * 10))
    # # Debugging: calculating total area under all curves:
    # totalArea = 0
    for i in range(nrows):
        for j in range(ncols):
            index = sectionID - 1 # The index for the distribution of section
            # Rescaling bin centres with box length
            if segParam.RESCALE_LENGTHS:
                binCentres_list[index] = binCentres_list[index] / boxLength
            # PlotSingleDistribution(axes[i][j], binCentres_list[index], distribution_list[index], f"Section {sectionID}")
            PlotSingleDistribution(axes[i][j], binCentres_list[index], number_density_list[index], f"Section {sectionID}")
            # # Debugging: Calculating area under the curve
            # area = plotMonomerDensity.CalculateAreaUnderCurve(distribution_list[index], binWidth)
            # totalArea += area
            # print(f"Section {sectionID} Area: {area} {axisLabel}")
            sectionID += 1
    # debugging:
    # print(f"Total area under all curves: {totalArea} {axisLabel}")
    if segParam.SHOW_TITLE:
        fig.suptitle("Monomer Distribution for %i %s polymer(s)\n with %i monomers each, %i sections, bin width = %.1lf" % (numberOfPolymers, architecture, numberOfMonomers, numberOfSections, binWidth))
    else: # Printing a shorter title
        fig.suptitle(segParam.GetPaperFiguresTitle(numberOfMonomers, architecture))
    xLabelModifier = ""
    if segParam.RESCALE_LENGTHS:
        xLabelModifier = r"/ $L$"
    fig.supxlabel(r"$z$ %s" % (xLabelModifier))
    if len(axisLabel) > 0:
        axisLabel = f"   [{axisLabel}]" # encasing in parenthesis if not empty
    fig.supylabel(r"$n(z$%s$)/(\pi R^2 dz)$%s" % (xLabelModifier, axisLabel))
    # fig.tight_layout()
    fig.subplots_adjust(wspace = 0.4)

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)
    # Saving figure:
    saveDirectory = "%sAnalysis/SectionalDistribution/" % (folder)
    Path(saveDirectory).mkdir(parents = True, exist_ok = True) # Creating the directory with parents
    saveFilePath = "%ssectional_monomer_distribution_r%i_s%i_b%02i%s" % (saveDirectory, runIndex, numberOfSections, round(binWidth*10), segParam.FIG_EXT)
    fig.savefig(saveFilePath)
    plt.close(fig)
    print(f"Plotted Section Monomer Distribution for {architecture} Run {runIndex}")

def PlotPolymerDistributions(showPlot: bool = False) -> None:
    """Plots the monomer distributions of each entire polymer in the whole cylinder."""
    fig, ax = plt.subplots(figsize = (12, 6))
    folder = GetFolder(numberOfMonomers, architecture, special_simulation=initializationProcedure)
    # Reading the regions monomer density: Adding them up to get the density for the entire polymer

    # Setting number of regions
    if reg.USE_REGIONS:
        numberOfRegions = reg.NUMBER_OF_REGIONS
    else:
        numberOfRegions = numberOfPolymers

    regionsPerPolymer = numberOfRegions // numberOfPolymers
    polymerDistDataframes = [] # List of summed distribution dataframes
    polymerProbDensities = [] # List of final prob densities for each polymer
    for polymerIndex in range(numberOfPolymers):
        regionDataframes = [] # List of distribution dataframes for the regions in the current polymer
        for i in range(regionsPerPolymer):
            regionIndex = polymerIndex * regionsPerPolymer + i
            filePath = GetDistributionFilePath(folder, runIndex, f"monomer_distribution_reg{regionIndex+1}.csv")
            # Reading the file:
            regionDataframes.append(pd.read_csv(filePath))
        # Adding distributions of regions to get polymer distribution:
        polymerDistDataframes.append(plotMonomerDensity.AddDistributions(regionDataframes))
        polymerProbDensities.append(list(polymerDistDataframes[polymerIndex].iloc[: , 1]))
    
    # Converting axial probability density into monomer number (volume) density:
    polymerNumberDensities = []
    for indexOfPolymer in range(numberOfPolymers):
        numberDensity = ConvertToMonomerNumberDensity(np.array(polymerProbDensities[indexOfPolymer]))
        polymerNumberDensities.append(numberDensity)

    # # Sanity check: Ensuring that the integral of each polymer density adds up to the number of monomers
    # binWidth = float(polymerDistDataframes[0].iloc[1, 0]) - float(polymerDistDataframes[0].iloc[0, 0])
    # for indexOfPolymer in range(numberOfPolymers):
    #     total_monomers = np.sum(polymerNumberDensities[indexOfPolymer]) * binWidth * (math.pi * diameter ** 2 / 4)
    #     print(f"Total monomers in Polymer {indexOfPolymer + 1}: {total_monomers} (Should be close to {numberOfMonomers})")
    #     # Adding up probability densities to check if they add up to 1
    #     total_probability = np.sum(polymerProbDensities[indexOfPolymer]) * binWidth
    #     print(f"Total probability in Polymer {indexOfPolymer + 1}: {total_probability} (Should be close to 1)")
    
    # Converting each of the probability distribution values to scientific notation
    polymerNumberDensities, yAxisLabel = pt.ConvertMultipleArraysToScientificNotation(polymerNumberDensities, (-1, 3), segParam.PREFERRED_ORDER)
    if len(yAxisLabel) > 0:
        yAxisLabel = f"   [{yAxisLabel}]"

    # Normalizing the bin Centres: this will make the binWidths different than what was chosen
    # TODO: Set number of Bins instead of binWidths? Won't be affected by scaling
    xLabelModifier = ""
    if segParam.RESCALE_LENGTHS:
        xLabelModifier = r"/ $L$"
    binWidth = float(polymerDistDataframes[0].iloc[1, 0]) - float(polymerDistDataframes[0].iloc[0, 0])
    # Debugging: Printing binWidth
    # print("Polymer Density binWidth:  %.12lf, binWidth label: %02i" % (binWidth, binWidth * 10))
    for indexOfPolymer in range(numberOfPolymers):
        # Plotting:
        binCentres = plotMonomerDensity.ShiftBinEdges(polymerDistDataframes[indexOfPolymer].iloc[:, 0])
        if segParam.RESCALE_LENGTHS:
            binCentres = binCentres / boxLength
        ax.plot(binCentres, polymerNumberDensities[indexOfPolymer], linestyle = '--', marker = '.', label = f"Polymer {indexOfPolymer + 1}")
        # Debugging: Calculating area under the curve
        # print(f"Polymer {indexOfPolymer+1} Area: {plotMonomerDensity.CalculateAreaUnderCurve(polymerDistributions[indexOfPolymer], binWidth)}")
    timeLabel = rf"{Decimal(numberOfSteps/TAU_0):.0E}$\tau_0$"
    if segParam.SHOW_TITLE:
        ax.set_title(f"Monomer Distribution for a {architecture} polymer\n Total: {numberOfPolymers} polymer(s), {numberOfMonomers} monomers each, {timeLabel}, Run {runIndex}")
    else: # Printing a shorter title
        ax.text(segParam.X_LIMS[0] * (1 - 0.1), segParam.Y_LIMS[1] * (1 - 0.1) , segParam.GetPaperFiguresTitle(numberOfMonomers, architecture), fontweight = 'bold', fontsize = 26, horizontalalignment = 'left')
        # ax.text(segParam.X_LIMS[0] * (1 - 0.1), segParam.Y_LIMS[1] * (1 - 0.1) , segParam.GetAliasSimulation(), fontweight = 'bold', fontsize = 26, horizontalalignment = 'left')
    ax.set_xlabel(r"$z$ %s" % (xLabelModifier), fontweight = 'bold')
    ax.set_ylabel(r"$n(z$%s$)/(\pi R^2 dz)$%s" % (xLabelModifier, yAxisLabel), fontweight = 'bold')
    leg = ax.legend(loc = 'lower center')

    # Setting major and minor ticks:
    # ax.xaxis.set_major_locator(MultipleLocator(0.2))
    # ax.xaxis.set_minor_locator(MultipleLocator(0.05))

    # Setting axis limits if necessary:
    if segParam.USE_COMMON_AXIS_LIMITS:
        ax.set_xlim(segParam.X_LIMS)
        ax.set_ylim(segParam.Y_LIMS)
    # Increasing legend line thickness:
    for line in leg.get_lines():
        line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())
    # Note: This monomer density is the probability of finding a monomer of a particular polymer 
    #       at a particular z position. In other words, summing up the densities of all the polymers
    #       gives the probability of finding a (any) monomer at a z position
    
    # Adjusting plot margins:
    if segParam.USE_MARGINS:
        fig.subplots_adjust(top = segParam.PLOT_MARGINS['top'], right = segParam.PLOT_MARGINS['right'], bottom = segParam.PLOT_MARGINS['bottom'])

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)
    # Saving figure:
    saveDirectory = "%sAnalysis/TotalDistribution/" % (folder)
    Path(saveDirectory).mkdir(parents = True, exist_ok = True) # Creating the directory with parents
    saveFilePath = "%stotal_monomer_distribution_r%i_b%02i%s" % (saveDirectory, runIndex, round(binWidth*10), segParam.FIG_EXT)
    fig.savefig(saveFilePath)
    plt.close(fig)
    print(f"Plotted Polymer Monomer Distribution for {architecture} Run {runIndex}")

def PlotRegionDistributions(showPlot: bool = False) -> None:
    """Plots the monomer distribution separately for each region in the entire cylinder"""
    fig, ax = plt.subplots()
    folder = GetFolder(numberOfMonomers, architecture, special_simulation=initializationProcedure)
    for regionIndex in range(reg.NUMBER_OF_REGIONS):
        filePath = GetDistributionFilePath(folder, runIndex, f"monomer_distribution_reg{regionIndex+1}.csv")
        binWidth = plotMonomerDensity.PlotDistribution(numberOfMonomers, architecture, filePath, reg.REGION_LABELS[regionIndex], ax, boxLength, runIndex = runIndex, numberOfSteps = numberOfSteps/TAU_0)
    # Note: This monomer density is the probability of finding a monomer of a particular region 
    #       at a particular z position. In other words, summing up the densities of all the regions
    #       gives the probability of finding a (any) monomer at a z position
    # Debugging: Printing binWidth
    # print("Regional Density binWidth:  %.12lf, binWidth label: %02i" % (binWidth, round(binWidth * 10)))
        
    # Adjusting plot margins:
    if segParam.USE_MARGINS:
        fig.subplots_adjust(top = segParam.PLOT_MARGINS['top'], right = segParam.PLOT_MARGINS['right'], bottom = segParam.PLOT_MARGINS['bottom'])

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)
    # Saving figure:
    saveDirectory = "%sAnalysis/RegionalDistribution/" % (folder)
    Path(saveDirectory).mkdir(parents = True, exist_ok = True) # Creating the directory with parents
    saveFilePath = "%sregion_monomer_distribution_r%i_b%02i%s" % (saveDirectory, runIndex, round(binWidth*10), segParam.FIG_EXT)
    fig.savefig(saveFilePath)
    plt.close(fig)
    print(f"Plotted Region wise Monomer Distribution for {architecture} Run {runIndex}")

def PlotSingleSnapshotRegionDistributions(showPlot: bool = False) -> None:
    """Plots the regional monomer distributions that have been computed from only a single snapshot"""
    # Plotting subplots: one for each polymer
    mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/subplots_bold.mplstyle")
    fig, ax = plt.subplots(nrows = numberOfPolymers, sharex = True, sharey = True)
    folder = GetFolder(numberOfMonomers, architecture, special_simulation=initializationProcedure)
    numberOfPolymerRegions = reg.NUMBER_OF_REGIONS // numberOfPolymers # Number of regions in one polymer
    for regionIndex in range(reg.NUMBER_OF_REGIONS):
        filePath = GetDistributionFilePath(folder, runIndex, f"monomer_distribution_reg{regionIndex+1}_single.csv")
        polymerIndex = regionIndex // numberOfPolymerRegions
        binWidth = plotMonomerDensity.PlotDistribution(numberOfMonomers, architecture, filePath, reg.REGION_LABELS[regionIndex], ax[polymerIndex], boxLength, runIndex = runIndex, numberOfSteps = 1/TAU_0, plotDensity = segParam.PLOT_SINGLE_SNAPSHOT_DENSITY)
        # Resetting title and axis labels:
        ax[polymerIndex].set_title("")
        ax[polymerIndex].set_xlabel("")
        ax[polymerIndex].set_ylabel("")
    # Setting title and axis labels:
    fig.suptitle(f"Monomer Distribution of 2 {architecture} polymers\n{numberOfMonomers} monomers each, 1 snapshot, Run {runIndex}")
    fig.supxlabel("Coordinate along cylinder axis (LJ Units)")
    if segParam.PLOT_SINGLE_SNAPSHOT_DENSITY:
        fig.supylabel("Probability Density")
    else:
        fig.supylabel("Number of Monomers")
    # Note: This monomer density is the probability of finding a monomer of a particular region 
    #       at a particular z position. In other words, summing up the densities of all the regions
    #       gives the probability of finding a (any) monomer at a z position
        
    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)
    # Saving figure:
    saveDirectory = "%sAnalysis/SingleSnapshotDistribution/" % (folder)
    Path(saveDirectory).mkdir(parents = True, exist_ok = True) # Creating the directory with parents
    fileName = single.GetFileNamePrefix()
    saveFilePath = "%s%s_b%02i_r%i%s" % (saveDirectory, fileName, int(round(binWidth*10)), runIndex, segParam.FIG_EXT)
    fig.savefig(saveFilePath)
    plt.close(fig)
    print(f"Plotted Single Snapshot Region wise Monomer Distribution for {architecture} Run {runIndex}")

def PlotAverageSingleSnapshotDistribution(showPlot: bool = False) -> None:
    """Plots the average probability distributions separately for each region
    The average is over a single snapshot for each of the runs"""
    folderPath = GetFolder(numberOfMonomers, architecture, special_simulation=initializationProcedure)
    totalProbList = []
    binsList= []
    binWidth = 0
    for regionIndex in range(reg.NUMBER_OF_REGIONS):
        for runIndex in range(1, numberOfRuns + 1):
            # reading file for a particular region
            filePath = GetDistributionFilePath(folderPath, runIndex, f"monomer_distribution_reg{regionIndex+1}_single.csv")
            df = pd.read_csv(filePath)
            binLeftEdges = np.array(df.iloc[: , 0])
            probDistribution = np.array(df.iloc[: ,1])
            # Adding probability distribution:
            if runIndex == 1: # first run; totalProbList not initialized for the region
                totalProbList.append(probDistribution) 
                # Adding the bins only once for each region:
                # calculating binWidth:
                if binWidth == 0:
                    binWidth = binLeftEdges[1] - binLeftEdges[0]
                binLastEdge = binLeftEdges[-1] + binWidth
                binEdges = np.append(binLeftEdges, binLastEdge)
                binsList.append(binEdges)

            else: # Add to existing numpy arrays in totalProbList
                totalProbList[regionIndex] += probDistribution
            # Adding binEdges:
    # Total probability should have been calculated
    # Normalizing with number of runs:
    for i in range(reg.NUMBER_OF_REGIONS):
        totalProbList[i] /= numberOfRuns

    # Plotting:
    title = f"Monomer Distribution of 2 {architecture} polymers\n{numberOfMonomers} monomers each, Averaged over {numberOfRuns} snapshots"
    fileName = ""
    fileName = single.GetFileNamePrefix()
    filePath = "%s/Analysis/SingleSnapshotDistribution/average_%s_b%02i%s" % (folderPath, fileName, int(round(binWidth*10)), segParam.FIG_EXT)
    single.PlotRegionalDistributions(binsList, totalProbList, runIndex, title, filePath, showPlot)
    print(f"Plotted average monomer density of single snapshots of {architecture} polymers over {numberOfRuns} runs.")

def PlotTotalDistribution(showPlot: bool = False) -> None:
    """Plots the monomer distribution of all polymers in the entire cylinder"""
    fig, ax = plt.subplots()
    folder = GetFolder(numberOfMonomers, architecture, special_simulation=initializationProcedure)
    distributions = []
    for regionIndex in range(1, reg.NUMBER_OF_REGIONS + 1):
        filePath = GetDistributionFilePath(folder, runIndex, f"monomer_distribution_reg{regionIndex}.csv")
        df = pd.read_csv(filePath)
        binLeftEdges = df.iloc[: , 0]
        binWidth = float(binLeftEdges[1]) - float(binLeftEdges[0])
        distribution = df.iloc[:, 1]
        distributions.append(np.array(distribution))
    binCentres = ShiftBinEdges(binLeftEdges, binWidth)
    
    # Adding the distributions:
    tot_distribution = np.zeros(len(distributions[0]))
    for distribution in distributions:
        tot_distribution += distribution
    
    # Converting axial probability density into monomer number (volume) density:
    numberDensity = ConvertToMonomerNumberDensity(np.array(tot_distribution))

    # Sanity check: Ensuring that the integral of the density adds up to the number of monomers
    total_monomers = np.sum(numberDensity) * binWidth * (math.pi * diameter ** 2 / 4)
    print(f"Total monomers: {total_monomers} (Should be close to {numberOfMonomers * numberOfPolymers})")
    # Adding up probability densities to check if they add up to 1
    total_probability = np.sum(tot_distribution) * binWidth
    print(f"Total probability: {total_probability} (Should be close to 1)")

    # Plotting:
    # Rescaling the x axis with box length
    xLabelModifier = ""
    if segParam.RESCALE_LENGTHS:
        binCentres = binCentres / boxLength
        xLabelModifier = r"/ $L$"
    # Converting to Scientific notation:
    numberDensity, axisLabel = pt.ConvertToScientificNotation(numberDensity, (-1, 3))
    ax.plot(binCentres, numberDensity, linestyle = '--', marker = '.')
    if segParam.SHOW_TITLE:
        ax.set_title(f"Total Monomer Distribution for a {architecture} polymer with {numberOfMonomers} monomers\n Total: {numberOfPolymers} polymer(s), {Decimal(numberOfSteps):.0E} timesteps, Run {runIndex}")
    else:
        ax.text(0, segParam.Y_LIMS[1] * (1 - 0.1), segParam.GetPaperFiguresTitle(numberOfMonomers, architecture), fontweight = 'bold', fontsize = 26, horizontalalignment = 'center')
    ax.set_xlabel(r"$z$ %s" % (xLabelModifier))
    if len(axisLabel) == 0:
        yLabelModifier = ""
    else:
        yLabelModifier = f" [{axisLabel}]"
    ax.set_ylabel(rf"$n(z/L)/(\pi R^2 dz)${yLabelModifier}")
    
    # Setting axis limits if necessary:
    if segParam.USE_COMMON_AXIS_LIMITS:
        ax.set_ylim(top = segParam.Y_LIMS[1])

    # Adjusting plot margins:
    if segParam.USE_MARGINS:
        fig.subplots_adjust(top = segParam.PLOT_MARGINS['top'], right = segParam.PLOT_MARGINS['right'], bottom = segParam.PLOT_MARGINS['bottom'])

    # Including the plot identification label if necessary:
    if segParam.INCLUDE_PLOT_IDENTIFICATION_LABEL:
        pt.IncludePlotIdentificationLabel(fig, segParam.ID_LABEL_TEXT)

    if showPlot:
        plt.show(block = True)

def PlotRadialDistribution(showPlot: bool = False) -> None:
    """Plots the radial distribution of the entire system in the cylinder and saves it as a figure"""
    folder = GetFolder(numberOfMonomers, architecture, special_simulation=initializationProcedure)
    filePath = GetDistributionFilePath(folder, runIndex, "radial_monomer_distribution.csv")
    radius = boxLength / segParam.ASPECT_RATIO / 2
    # Plotting:
    fig, ax = plt.subplots()

    # Marking the section radius threshold:
    radiusThreshold = (radius - 0.5) / math.sqrt(2)
    if segParam.RESCALE_RADIUS:
        radiusThreshold = radiusThreshold / radius
    ax.axvline(x = radiusThreshold, linestyle = '--', color = 'r', label = "Section threshold")

    binWidth = plotMonomerDensity.PlotRadialDistribution(numberOfMonomers, architecture, filePath, "", ax, radius, runIndex, numberOfSteps = numberOfSteps / TAU_0)
    leg = ax.legend(loc = 'lower left')

    # Setting axis limits if necessary:
    if segParam.USE_COMMON_AXIS_LIMITS:
        ax.set_ylim(top = segParam.Y_LIMS[1])

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

    # Saving figure:
    saveDirectory = "%sAnalysis/RadialDistribution/" % (folder)
    Path(saveDirectory).mkdir(parents = True, exist_ok = True) # Creating the directory with parents
    saveFilePath = "%sradial_monomer_distribution_r%i_b%02i%s" % (saveDirectory, runIndex, int(round(binWidth*10)), segParam.FIG_EXT)
    fig.savefig(saveFilePath)
    plt.close(fig)
    print(f"Plotted Radial Monomer Distribution for {architecture} Run {runIndex}")

def PlotAndSaveRegionwiseRadialDistribution(showPlot: bool = False) -> None:
    """Plots the radial monomer distributions of all the regions in the system and saves the figure as an image."""
    fig, ax = plt.subplots()
    folder = GetFolder(numberOfMonomers, architecture, special_simulation=initializationProcedure)
    for regionIndex in range(1, reg.NUMBER_OF_REGIONS + 1):
        # getting data from CSV file:
        filePath = f"{folder}run{runIndex}/monomer_distribution/radial_monomer_distribution_reg{regionIndex}.csv"
        radius = boxLength / segParam.ASPECT_RATIO / 2
        binWidth = plotMonomerDensity.PlotRadialDistribution(numberOfMonomers, architecture, filePath, reg.REGION_LABELS[regionIndex-1], ax, radius, runIndex)
    leg = ax.legend()
    # Increasing legend line thickness:
    for line in leg.get_lines():
        line.set_linewidth(pt.LEGEND_LINEWIDTH * line.get_linewidth())
    if showPlot:
        plt.show(block = True)
    fig.savefig(f"{folder}Analysis/RadialDistribution/regionwiseDistribution_r{runIndex}.png")
    plt.close(fig) # closing the figure to save memory
    print(f"Plotted regionwise Radial Monomer Distribution for {architecture} Run {runIndex}")
 
#script:
if __name__ == "__main__":
    SetArguments()
    if runIndex == -1:
        for runIndex in range(1, numberOfRuns + 1):
            mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/subplots_bold.mplstyle") # Setting the mpl style sheet
            PlotDistributionForAllSections(numberOfSections = 8, showPlot = False)
            mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/bold.mplstyle") # Setting the mpl style sheet
            if reg.USE_REGIONS:
                if not segParam.USE_SINGLE_SNAPSHOT:
                    PlotRegionDistributions(False)
                else:
                    PlotSingleSnapshotRegionDistributions(False)
            if segParam.PLOT_INIT_COMPARISON:
                plotMonomerDensity.PlotRadialInitializationComparison(segParam.SPECIAL_SIMULATIONS, numberOfMonomers, architecture, boxLength, runIndex, sysPaths.CREATE_INITIAL_STATES, False)
            else:
                PlotRadialDistribution(False)
            PlotPolymerDistributions(False)

        print(f"Successfully plotted the Monomer Distributions for the {numberOfRuns} runs of {architecture}")
    elif runIndex == 0:
        if not segParam.USE_SINGLE_SNAPSHOT:
            print("The passed run index 0 is meaningless unless using the user wants to calculate the distribution from a single snapshot.")
        else:
            PlotAverageSingleSnapshotDistribution(True)
    else:
        mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/subplots_big_bold.mplstyle") # Setting the mpl style sheet
        # PlotDistribution(True)
        PlotDistributionForAllSections(numberOfSections = 8, showPlot = True)
        mpl.style.use(f"{sysPaths.GLOBAL_SCRIPTS}Plotting_Styles/big_bold.mplstyle") # Setting the mpl style sheet
        if reg.USE_REGIONS:
            if not segParam.USE_SINGLE_SNAPSHOT:
                PlotRegionDistributions(True)
            else:
                PlotSingleSnapshotRegionDistributions(True)
        if segParam.PLOT_INIT_COMPARISON:
            plotMonomerDensity.PlotRadialInitializationComparison(segParam.SPECIAL_SIMULATIONS, numberOfMonomers, architecture, boxLength, runIndex, sysPaths.CREATE_INITIAL_STATES, True)
        else:
            PlotRadialDistribution(True)
        # PlotAndSaveRegionwiseRadialDistribution(True)
        PlotPolymerDistributions(True)
        # PlotTotalDistribution(True)
