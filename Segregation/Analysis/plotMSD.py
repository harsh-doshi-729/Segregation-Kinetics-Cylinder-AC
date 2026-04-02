# Plots the mean squared displacement (MSD) calculated from the simulation data

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import Segregation_Parameters as segParam
from pathlib import Path

from typing import Tuple, List
from numpy.typing import ArrayLike
# importing the scripts that contains all the local paths
# The path to the directory containing this module is set as the environment variable PYTHONPATH
import system_file_paths as sysPaths
import PlottingTools as pt

# Global variables:
numberOfPolymers = 2
numberOfMonomers = 200
architecture = ""
runIndex = -1

MSD_mode = "" # The mode indicating how the MSD has been calculated
acceptedModes = ["COM", "relative_COM", "z_COM"]

def SetConstants() -> None:
    """
    Reading arguments from the command line and setting the global variables accordingly.
    """
    global numberOfMonomers
    global architecture
    global runIndex
    global MSD_mode

    numberOfMandatoryArguments = 3 # numberOfMonomers, architecture, mode
    if len(sys.argv) < numberOfMandatoryArguments + 1: # +1 for the script name
        print(f"Not enough arguments provided. Please pass the number of monomers, the architecture, and the MSD mode as command line arguments.")
        print("Usage: python plotMSD.py <numberOfMonomers> <architecture> <mode> [runIndex]")
        print("An optional argument for the run index may be passed after the mandatory arguments.")
        print(f"Accepted MSD modes: {acceptedModes}")
        print("Example: python plotMSD.py 200 Arc2 COM 1")
        sys.exit(1)
    else:
        try:
            numberOfMonomers = int(sys.argv[1])
        except ValueError:
            print("The number of monomers must be an integer. Terminating.")
            sys.exit(1)
        architecture = sys.argv[2]
        MSD_mode = sys.argv[3]
        if MSD_mode not in acceptedModes:
            print(f"The passed MSD mode '{MSD_mode}' is not recognized. Accepted modes are: {acceptedModes}. Terminating.")
            sys.exit(1)
        # Optional arguments:
        if len(sys.argv) > numberOfMandatoryArguments + 1:
            try:
                runIndex = int(sys.argv[numberOfMandatoryArguments + 1])
            except ValueError:
                print("The run index must be an integer. Terminating.")
                sys.exit(1)

def GetArchFolder() -> str:
    """
    Returns the path to the folder containing the data for the given architecture.
    The path is constructed based on the global variables set from command line arguments.
    """
    baseFolder = sysPaths.GetFolder(sysPaths.NEW_SEGREGATION, numberOfMonomers, sysPaths.SPECIAL_SIMULATION)
    archFolder = f"{baseFolder}{architecture}/"
    if not Path(archFolder).is_dir():
        print(f"The architecture folder '{archFolder}' does not exist. Please check the parameters and try again. Terminating.")
        sys.exit(1)
    return archFolder

def GetDataFilePath(runIndex: int) -> str:
    """
    Returns the path to the data file containing the MSD data.
    The path is constructed based on the global variables set from command line arguments.
    """
    dataFilePath = f"{GetArchFolder()}run{runIndex}/{MSD_mode}_MSD.csv"
    return dataFilePath

def ReadMSD(filePath: str) -> Tuple[ArrayLike, ArrayLike]:
    """
    Reads the MSD data from the file and returns the time and MSD arrays.
    The data file is expected to be in CSV format with two columns: (delta)time and MSD.
    Args:
        -filePath: The path to the data file.
    Returns:
        -time: The array of time values.
        -msd: The array of MSD values.
    """
    if not Path(filePath).is_file():
        print(f"The data file '{filePath}' does not exist. Please check the path and try again. Terminating.")
        sys.exit(1)
    
    # Reading the data
    try:
        data = pd.read_csv(filePath)
        time = data.iloc[:, 0].to_numpy()
        msd = data.iloc[:, 1].to_numpy()
        return time, msd
    except Exception as e:
        print(f"An error occurred while reading the MSD data file: {e}")
        sys.exit(1)

# def ReadAndPlotMSD(ax: mpl.axes.Axes, MSD_mode: str, runIndex: int) -> None:
#     """
#     Reads the MSD data from the file and plots it on the given axes.
#     The data file is expected to be in CSV format with two columns: (delta)time and MSD.
#     Args:
#         -ax: The matplotlib axes on which to plot the data.
#         -MSD_mode: The mode indicating how the MSD has been calculated.
#         -runIndex: The index of the run to read the data from.
#     """
#     dataFilePath = GetDataFilePath(runIndex)
#     time, msd = ReadMSD(dataFilePath)
    
#     # Plotting the data
#     ax.plot(time, msd)

def CalculateAverageMSD(MSD_mode: str) -> Tuple[ArrayLike, ArrayLike]:
    """
    Calculates the average MSD over all runs.
    Args:
        -MSD_mode: The mode indicating how the MSD has been calculated.
    Returns:
        -time: The array of time values.
        -avg_msd: The array of average MSD values.
    """ 
    all_msds = []
    time = None
    for run in range(1, segParam.NUMBER_OF_RUNS + 1):
        dataFilePath = GetDataFilePath(run)
        t, msd = ReadMSD(dataFilePath)
        if time is None:
            time = t
        all_msds.append(msd)
    
    # Calculating the average MSD
    avg_msd = np.mean(np.array(all_msds), axis=0)
    return time, avg_msd

def FitLinearMSD(time: ArrayLike, msd: ArrayLike) -> Tuple[float, float]:
    """
    Fits a straight line (y = a*x + b) to the MSD data over the given time range.
    Returns:
        - slope (a): The fitted slope.
        - intercept (b): The fitted intercept.
    """
    coeffs = np.polyfit(time, msd, 1)
    slope, intercept = coeffs
    return slope, intercept

def FitPowerLawMSD(time: ArrayLike, msd: ArrayLike) -> Tuple[float, float]:
    """
    Fits a power law (y = a * x^b) to the MSD data.
    Returns:
        - exponent (b): The fitted exponent.
        - prefactor (a): The fitted prefactor.
    """
    # Avoid log(0) by filtering out zero or negative values
    mask = (time > 0) & (msd > 0)
    log_time = np.log(time[mask])
    log_msd = np.log(msd[mask])
    coeffs = np.polyfit(log_time, log_msd, 1)
    exponent = coeffs[0]
    log_prefactor = coeffs[1]
    prefactor = np.exp(log_prefactor)
    return exponent, prefactor

def Plot_MSD(showPlot: bool = False) -> None:
    """
    Plots the center of mass (COM) MSD for the given architecture and number of monomers.
    Args:
        -showPlot: A flag to indicate whether the plot should be launched in the interactive mode
    """
    mpl.style.use(f"{sysPaths.POLYMER_PHYSICS}Scripts/Plotting_Styles/bold.mplstyle")
    if segParam.USE_LOGLOG_PLOT:
        mpl.rcParams['axes.formatter.limits'] = (-3, 8)
    fig, ax = plt.subplots()

    if segParam.USE_LOGLOG_PLOT:
        plottingFunction = ax.loglog
    else:
        plottingFunction = ax.plot

    if runIndex == -1: # Average over all runs
        time, msd = CalculateAverageMSD(MSD_mode)
        runLabel = f"Averaged over {segParam.NUMBER_OF_RUNS} runs"
    else: # Specific run
        time, msd = ReadMSD(GetDataFilePath(runIndex))
        runLabel = f"Run {runIndex}"
    time = time / segParam.TAU_0    
    time, axisLabel = pt.ConvertToScientificNotation(time)
    plottingFunction(time, msd, marker = '.')
    # Fitting power law:
    if sysPaths.IsCylinderInfinite():
        exponent, prefactor = FitPowerLawMSD(time, msd)
    else:
        exponent, prefactor = (1, 0.01)
    fit_msd = prefactor * time ** exponent
    ax.plot(time, fit_msd, linestyle='--', color = 'r', label = f"Exponent = {exponent:.2f}")
    ax.legend()

    # Setting the plot labels and title
    ax.set_xlabel(rf"$\Delta$Time ({axisLabel}$\tau_0$)")
    ax.set_ylabel(rf"MSD ($\sigma^2$)")
    ax.set_title(f"{MSD_mode} MSD for {numberOfPolymers} {architecture} polymers\n{numberOfMonomers} monomers each, {runLabel}")

    if showPlot:
        plt.show(block = True)
    
    # Saving figure:
    saveFolder = f"{GetArchFolder()}Analysis/MSD/"
    Path(saveFolder).mkdir(parents=True, exist_ok=True)
    if runIndex == -1:
        fileRunLabel = "all"
    else:
        fileRunLabel = f"r{runIndex}"
    fig.savefig(f"{saveFolder}{MSD_mode}_MSD_{fileRunLabel}{segParam.FIG_EXT}")
    plt.close(fig)
    # Printing success message
    print(f"The {MSD_mode} MSD plot {runLabel} was plotted and saved successfully.")

if __name__ == "__main__":
    SetConstants()
    if MSD_mode == "COM":
        Plot_MSD(showPlot=True)
    elif MSD_mode == "relative_COM":
        Plot_MSD(showPlot=True)
    elif MSD_mode == "z_COM":
        Plot_MSD(showPlot=True)