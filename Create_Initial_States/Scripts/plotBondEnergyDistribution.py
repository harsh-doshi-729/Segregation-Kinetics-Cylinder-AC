import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import sys
import math

# importing the scripts that contains all the local paths
# The path to the directory containing this module is set as the environment variable PYTHONPATH
import system_file_paths as sysPaths 
import Segregation_Parameters as segParam

numberOfMonomers = 200
numberOfPolymers = 2
architecture = "Arc0"
runIndex = -1
numberOfRuns = 50

numberOfSections = 8
directory = sysPaths.CREATE_INITIAL_STATES
if sysPaths.USE_R_G_FOLDER:
    directory = sysPaths.R_G


def SetConstants():
    """Reads the command line arguments and sets constants"""
    global numberOfMonomers
    global architecture
    global runIndex
    global numberOfPolymers
    
    numberOfMandatoryArguments = 2
    if len(sys.argv) <= numberOfMandatoryArguments:
        print(f"Not enough arguments passed! Please pass the number of monomers and the architecture as arguments to the script.")
        print("For example: script.py 200 Arc0")
        print("An optional argument for a run number can be passed to plot for just that run.")
        quit(1)
    
    try:
        numberOfMonomers = int(sys.argv[1])
        architecture = sys.argv[2]

        if len(sys.argv) > numberOfMandatoryArguments + 1:
            runIndex = int(sys.argv[numberOfMandatoryArguments+1])

        if sysPaths.USE_R_G_FOLDER:
            numberOfPolymers = 1
    except ValueError:
        print("Either the number of monomers or the run number could not be converted to an integer. Terminating")



def PlotDistribution(binCentres: list[float], distribution: list[float], ax, title: str = "") -> None:
    """Plots the distribution in the passed pyplot axis"""
    ax.plot(binCentres, distribution, linestyle = '--', marker = '.')
    ax.set_yscale('log')
    ax.tick_params(axis = 'both', labelleft = True, labelbottom = True) # ensuring that the usual x and y tick labels are visible
    if len(title) > 0:
        ax.set_title(title)
    
def ShiftBinEdges(binEdges):
    """Returns a numpy array of the N bin centres based on the N+1 bin edges passed"""
    binWidth = binEdges[1] - binEdges[0]
    binCentres = binEdges[:-1] + 0.5 * binWidth
    return binCentres

def GetDistributionFolder(baseFolder: str, numberOfMonomers: int, architecture: str, runIndex: int) -> str:
    """Returns the folder path containing the bond energy files by taking special_simulation into account"""
    if not sysPaths.USE_R_G_FOLDER:
        folder = sysPaths.GetFolder(baseFolder, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)
        suffix = f"{architecture}/run{runIndex}/bond_energy/"
    else:
        folder = sysPaths.GetRgFolder(baseFolder, numberOfMonomers, sysPaths.R_G_SIMULATION)
        suffix = f"{architecture}/bond_energy/"
    return f"{folder}{suffix}"

def ReadAndPlotDistributions(numberOfSections: int, showPlot: bool = False):
    """Reads the various distribution files and plots the distribution in a collage of plots"""
    binCentresList = []
    distributionsList = []
    numberOfDataPoints = []

    energy_range = (0, 2) # in terms of kBT
    bins = 100

    folderPath = GetDistributionFolder(directory, numberOfMonomers, architecture, runIndex)
    for i in range(1, numberOfSections+1):
        filePath = f"{folderPath}bond_energies_sect{i}.csv"
        df = pd.read_csv(filePath)
        # Third column contains the bond energies
        bondEnergies = df.iloc[: , 2]
        numberOfDataPoints.append(len(bondEnergies))
        distribution, binEdges = np.histogram(bondEnergies, bins = bins, range = energy_range, density = False)

        # debugging:
        # if i == 1:
            # print(f"Distribution: {distribution}")

        # Adding histogram:
        binCentresList.append(ShiftBinEdges(binEdges))
        distributionsList.append(distribution)
    # Normalizing the distribution to get probability across single sections:
    totalDataPoints = sum(numberOfDataPoints)
    # print(f"Total Number of Datapoints = {totalDataPoints}")
    for i in range(len(distributionsList)):
        distributionsList[i] = distributionsList[i] / numberOfDataPoints[i]

    # debugging:
    # print(f"Distribution: {distributionsList[0]}")

    # plotting the distributions:
    nrows = 2
    ncols = numberOfSections // nrows
    fig, axes = plt.subplots(nrows = nrows, ncols = ncols, sharey = True, sharex = True, figsize = (20, 10))
    sectionID = 0
    for i in range(nrows):
        for j in range(ncols):
            PlotDistribution(binCentresList[sectionID], distributionsList[sectionID], axes[i][j], f"Section {sectionID+1} (N = {numberOfDataPoints[sectionID]})")
            sectionID += 1

    fig.suptitle(f"Bond Energy Distribution for {numberOfPolymers} {architecture} polymer(s)\n{numberOfMonomers} monomers each, {numberOfSections} sections")
    fig.supylabel("Probability")
    fig.supxlabel("Bond Energy (kBT)")
    # fig.tight_layout()
    fig.subplots_adjust(wspace = 0.3)

    if showPlot:
        plt.show(block = True)
    # saving figure:
    fig.savefig(f"{folderPath}bond_energy_distribution_s{numberOfSections}.png")

def PlotBondEnergies(ax, bondEnergies):
    """Plots the distribution of the passed bond energies and fits the ideally expected curve"""
    # mean of bond energies:
    mean = np.mean(bondEnergies)
    print("Mean Bond Energy = %.3lf" % (mean))

    distribution, binEdges = np.histogram(bondEnergies, bins = 100, density = True)
    binCentres = ShiftBinEdges(binEdges)
    # Plotting
    ax.plot(binCentres, distribution, linestyle = '', marker = '.')
    ax.set_title(f"Bond Energy Distribution for {numberOfPolymers} {architecture} polymer(s)\n{numberOfMonomers} monomers each, N = {len(bondEnergies)}")
    ax.set_xlabel("Bond Energy (kBT)")
    ax.set_ylabel("Probability Density")
    ax.set_yscale('log')


    # fitting the modified exponential (as expected from tranforming the bond length random variable):
    # ideal distribution density: sqrt(K / PI) / sqrt(x) * exp(-x); note the non linear scaling prefactor
    # endIndex = np.where(distribution == 0)[0][0]
    # endIndex = np.where(binCentres >= 2)[0][0] # taking data only till 2 kBT
    # logDistribution = np.log(distribution[:endIndex])
    # startIndex = 10
    # slope, intercept = np.polyfit(binCentres[startIndex:endIndex], logDistribution[startIndex:], deg = 1)
    # # slope is k, the bond constant
    # print("Fitting parameters: Slope = %.2lf, Intercept = %.6lf" % (slope, intercept))
    # plotting with ideal curve:
    slope = -1
    intercept = 0
    startIndex = 0
    endIndex = -1
    # predicted_y = np.exp(slope * binCentres) * math.exp(intercept)
    bond_coeff = 1000
    predicted_y = 1 / math.sqrt(math.pi) / np.sqrt(binCentres) * np.exp(-binCentres)
    ax.plot(binCentres[startIndex:endIndex], predicted_y[startIndex: endIndex], linestyle = '--', color = 'r', label = "Ideal expected curve")
    # plt.plot(binCentres[startIndex:endIndex], predicted_y[startIndex: endIndex], linestyle = '--', color = 'r', label = "Boltzmann distribution\nSlope = %.2lf" % (slope))
    ax.legend()

def ReadAndPlotBondEnergyDistribution(showPlot: bool = False) -> None:
    """Reads the bond energies for the entire cylinder and plots it as a single distribution"""
    folder = sysPaths.GetFolder(directory, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)
    bondFolderPath = GetDistributionFolder(directory, numberOfMonomers, architecture, runIndex)
    bondEnergies = []
    
    filePath = f"{bondFolderPath}bond_energies.csv"
    df = pd.read_csv(filePath)
    # Third column contains the bond energies
    bondEnergies += list(df.iloc[: , 2])
    
    fig, ax = plt.subplots()
    PlotBondEnergies(ax, bondEnergies)

    if showPlot:
        plt.show(block = True)

    fig.savefig(f"{folder}{architecture}/Analysis/bond_energy/bond_energy_distribution_r{runIndex}.png")
    plt.close(fig)
    print(f"Plotted Bond Energy Distribution for {architecture} Run {runIndex}")


def ReadAndPlotSystemEnergy(showPlot: bool = False) -> None:
    """Reads the thermo data file generated by a LAMMPS simulation and plots the enery distribution"""
    Rg_directory = sysPaths.R_G
    folderPath = f"{sysPaths.GetFolder(Rg_directory, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)}{architecture}/"

    filePath = f"{folderPath}thermo_output.dat"
    df = pd.read_csv(filePath, sep = ' ')
    ke = df.iloc[:, 1]
    pe = df.iloc[:, 2]
    totEnergy = df.iloc[:, 3]

    data = [ke, pe, totEnergy]
    binCentresList = []
    distributionsList = []
    for i in range(3):
        distribution, binEdges = np.histogram(data[i], bins = 'fd', density = True)
        distributionsList.append(distribution)
        binCentresList.append(ShiftBinEdges(binEdges))
    
    # Plotting:
    fig, ax = plt.subplots()
    labels = ["KE", "PE", "Total Energy"]
    for i in range(3):
        ax.plot(binCentresList[i], distributionsList[i], linestyle = '--', marker = '.', label = labels[i])

    # Fitting gaussian for total energy:
    mean = np.mean(data[2])
    std = np.std(data[2])
    y = GaussianPDF(binCentresList[2], mean, std)
    ax.plot(binCentresList[2], y, linestyle = '--', color = 'red', label = "Gaussain Fit\nMean = %.4lf\nSTD = %.4lf\n" % (mean, std))
    ax.legend()
    ax.set_title(f"System Energy Distribution for {numberOfPolymers} {architecture} polymer(s)\n{numberOfMonomers} monomers each")
    ax.set_ylabel("Probability Density")
    ax.set_xlabel("System Energy (in kBT units)")

    if showPlot:
        plt.show(block = True)




def GaussianPDF(x, mean: float = 0, std: float = 1):
    """Returns the value(s) of the gaussian distribution (Probability density function) for a passed value of the quantity"""
    return 1 / (math.sqrt(2 * math.pi) * std) * np.exp(-(x - mean)**2 / (2 * std**2))

def PlotBondLengths(ax, bondLengths):
    # mean of bond lengths:
    mean = np.mean(bondLengths)
    print("Mean Bond Length = %.4lf" % (mean))
    # standard deviation:
    std = np.std(bondLengths)
    print("Standard deviation = %.4lf" % (std))

    distribution, binEdges = np.histogram(bondLengths, bins = 'fd', density = True) 
    binCentres = ShiftBinEdges(binEdges)
    # Plotting
    ax.plot(binCentres, distribution, linestyle = '', marker = '.', label = "Data")
    ax.set_title(f"Bond Length Distribution for {numberOfPolymers} {architecture} polymer(s)\n{numberOfMonomers} monomers each, N = {len(bondLengths)}")
    ax.set_xlabel("Bond Length (simulation units)")
    ax.set_ylabel("Probability Density")

    # fitting Gaussian:
    y = GaussianPDF(binCentres, mean, std)
    ax.plot(binCentres, y, linestyle = '--', color = 'red', label = "Gaussian fit\nMean = %.4lf\nStd = %.4lf" % (mean, std))
    ax.legend()

def ReadAndPlotBondLengths(showPlot: bool = False) -> None:
    """Reads the bond lengths for the entire cylinder and plots it as a single distribution"""
    folderPath = GetDistributionFolder(directory, numberOfMonomers, architecture, runIndex)
    bondLengths = []
    
    filePath = f"{folderPath}bond_lengths.csv"
    df = pd.read_csv(filePath)
    # Third column contains the bond lengths
    bondLengths += list(df.iloc[: , 2])
    
    fig, ax = plt.subplots()
    PlotBondLengths(ax, bondLengths)

    if showPlot:
        plt.show(block = True)
    
    fig.savefig(f"{folderPath}bond_length_distribution.png")
    plt.close(fig)
    print(f"Plotted Bond Lengths Distribution for {architecture} Run {runIndex}")

if __name__ == "__main__":
    SetConstants()
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}/Scripts/Plotting_Styles/bold.mplstyle")
    if runIndex == -1 and not sysPaths.USE_R_G_FOLDER:
        for runIndex in range(1, numberOfRuns+1):
            # ReadAndPlotDistributions(numberOfSections = 10, showPlot = True)
            ReadAndPlotBondEnergyDistribution(showPlot = False)
            # ReadAndPlotBondLengths(True)
            # ReadAndPlotSystemEnergy(True)
    else:
        # ReadAndPlotBondEnergyDistribution(showPlot = True)
        ReadAndPlotBondLengths(showPlot = True)