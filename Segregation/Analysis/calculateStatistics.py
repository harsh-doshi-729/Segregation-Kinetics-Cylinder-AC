import pandas as pd
import numpy as np
import sys
import math
import system_file_paths as sysPaths 
import Segregation_Parameters as segParam
import plotSegTimeDistribution

from typing import List, Tuple

# polymer info:
architecture = "Arc2"
numberOfPolymers = 2
numberOfMonomers = 200

# file info:
directory = f"{sysPaths.POLYMER_PHYSICS}LAMMPS_runs/cluster-data/new_segregation/" # without any special simulation folder
special_simulation = ""

def SetConstants():
    """Reads the argument(s) passed while invoking the script from the command line and sets it to the architecture and numberOfMonomers"""
    global architecture
    global numberOfMonomers
    global special_simulation
    
    numberOfMandatoryArguments = 2
    if len(sys.argv) < numberOfMandatoryArguments:
        print("Not enough arguments specified while invoking the script! Please pass the name of the architecture and number of monomers in the following format:")
        print("python /<path>/plotCoMDistribution.py <noOfMonomers> <architecture>")
        print("An optional argument for the special simulation can be passed. No need to set its value in the config files")
        quit() # terminating the script
    
    architecture = sys.argv[2]
    numberOfMonomers = sys.argv[1]
    # checking if the number of monomers is valid:
    if not numberOfMonomers.isdigit():
        print(f"The entered numberOfMonomers {numberOfMonomers} cannot be converted to a number! Please provide a valid number")
        quit()

    numberOfMonomers = int(numberOfMonomers)

    if len(sys.argv) > numberOfMandatoryArguments + 1: # optional argument passed:
        special_simulation = sys.argv[numberOfMandatoryArguments+1]

def CalculateWithoutOutliers(segTimes: List[float]) -> Tuple[float, float, List[float]]:
    """Calculates the mean and standard deviation of the passed times after excluding outliers.
    Returns the mean, std. dev., and a list of detected outliers"""

    # Algorithm for detecting outliers:
    # calculate interquartile range (IQR)
    # define upper and lower "fences" as Q1-(1.5*IQR) and Q3+(1.5*IQR)
    # A point is considered an outlier if it lies beyond these fences
    segTimes.sort()
    Q1 = np.percentile(segTimes, 25) # First Quartile
    Q3 = np.percentile(segTimes, 75) # Third quartile
    IQR = Q3 - Q1
    factor = 2
    fences = [Q1 - factor * IQR, Q3 + factor * IQR]
    mean = 0
    std = 0
    outliers = []
    for segTime in segTimes:
        if segTime > fences[0] and segTime < fences[1]: # not an outlier
            mean += segTime
        else:
            outliers.append(segTime)

    numberOfDataPoints = len(segTimes) - len(outliers) # excluding data points
    mean = mean / numberOfDataPoints
    for segTime in segTimes:
        if segTime > fences[0] and segTime < fences[1]: # not an outlier
            std += (segTime - mean) ** 2
    
    std = math.sqrt(std / (numberOfDataPoints - 1))
    print(f"Q1 = {Q1}; Q3 = {Q3}; IQR = {IQR}")
    return mean, std, outliers

def GetFolder() -> str:
    """Returns the folder path where the segregation times csv file is located"""
    global special_simulation
    if len(special_simulation) != 0:
        special_simulation = f"{special_simulation}/"
    return f"{directory}{special_simulation}b{numberOfMonomers}/{architecture}/Analysis/"

def CalculateAndPrintStatistics() -> None:
    print(f"Printing segregation time statistics for {architecture};", end = " ")
    if len(special_simulation ) != 0:
        print(f"Special Simulation: {special_simulation}")
    else:
        print("")
    

    segTimes, segregation_criterion = plotSegTimeDistribution.ReadSegTimes(GetFolder())
    segTimes = plotSegTimeDistribution.RescaleTimesDict(segTimes, segParam.TAU_0)
    print(f"Segregation criterion: {segregation_criterion}")
    mean, std, outliers = CalculateWithoutOutliers(list(segTimes.values()))
    print(f"Mean = {mean}\n Std. Dev. = {std}")
    print(f"Excluded Outliers: {outliers}")

if __name__ == "__main__":
    SetConstants()
    CalculateAndPrintStatistics()
