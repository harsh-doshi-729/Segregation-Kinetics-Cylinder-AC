import numpy as np
import sys
import matplotlib.pyplot as plt

# importing the scripts that contains all the local paths
# The path to the directory containing this module is set as the environment variable PYTHONPATH
import system_file_paths as sysPaths 
import Segregation_Parameters as segParam

# This script reads the Average CoM Distances calculated for all the mixed states for different architectures, and plots them against the architecture
directory = sysPaths.CREATE_INITIAL_STATES
architectures = segParam.AOI_LIST
numberOfMonomers = 0
box_length = 25

def SetConstants():
    """Reads the arguments passed to the script and assigns their values to corresponding variables"""
    global numberOfMonomers
    global directory
    global box_length
    if len(sys.argv) < 2:
        print("Not enough arguments passed. Please pass the number of monomers from the command line.")
        quit()
    else:
        numberOfMonomers = int(sys.argv[1])
        directory = f"{directory}b{numberOfMonomers}/"
        if numberOfMonomers == 500:
            box_length = 35

def ReadAverageCoMDistance(architecture: str) -> float:
    """Reads the Average CoM Distance file in the specified architecture subdirectory and returns the read CoM Distance."""
    fileName = "AverageCoMDistances.log"
    filePath = f"{directory}{architecture}/{fileName}"
    file = open(filePath, "r")
    lines = file.readlines()
    # extracting the average CoM distance over all runs:
    averageCoMDistance = float(lines[-1].split(" ")[-1])
    print(f"Average CoM Distance for {architecture}: {averageCoMDistance}")
    return averageCoMDistance

def GetDistanceToSegregation(architecture: str) -> float:
    """Reads the Average CoM Distance from the particular architecture file; calculates and returns the distance the CoM Distance must increase by to be called segregated.
    This quantity serves as a parameter to signify how far away the mixed state is from segregation."""
    averageCoMDistance = ReadAverageCoMDistance(architecture)
    distanceToSegregation = box_length / 2 - averageCoMDistance
    return distanceToSegregation


def PlotAverageCoMDistances():
    """Reads and plots the average CoM Distance for all the architectures"""
    averageCoMDistances = []
    for architecture in architectures:
        averageCoMDistances.append(ReadAverageCoMDistance(architecture))
    
    # calculating another parameter: a measure of the distance that the mixed state must evolve to be called segregated:
    distancesToSegregation = box_length / 2 - np.array(averageCoMDistances)

    fig, ax = plt.subplots()
    # ax.scatter(architectures, averageCoMDistances)
    ax.scatter(architectures, distancesToSegregation)
    ax.set_title("Average Distance that must be crossed by polymers to segregate\nAveraged over 10⁷ iterations and 50 runs")
    # ax.set_title("Average CoM Distance between both polymers in mixed state\nAveraged over 10⁷ iterations and 50 runs")
    ax.set_xlabel("Architecture")
    ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation = 45)
    # ax.set_ylabel("Average CoM Distance")
    ax.set_ylabel("Box Length / 2 - Average CoM Distance")
    # save figure:
    fig.savefig(f"{directory}DistanceToSegregation.png")

# Script:
if __name__ == "main":
    SetConstants()
    PlotAverageCoMDistances()
