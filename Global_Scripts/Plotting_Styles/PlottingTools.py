# This is a helper script that contains a collection of tools useful for plotting with matplotlib.pyplot

import numpy as np
import pandas as pd
import matplotlib as mpl
import math
import re # Regular expressions
from numpy.typing import ArrayLike
from typing import Tuple, List
import sys
import os

import system_file_paths as sysPaths
import Segregation_Parameters as segParam

LEGEND_LINEWIDTH = 1.5 # A factor to be multiplied to the line width of the legend while plotting
# Goal is to make the legend lines more visible (thicker)

def ConvertToScientificNotation(array: ArrayLike, limits: None|tuple[int, int] = None, preferredOrder: None|int = None) -> Tuple[ArrayLike, str]:
    """Accepts a list/iterable of value and returns a tuple. 

    Args:
        array: A 1D array containing the numerical data to be converted
        limits: Pair (tuple) of integers to set lower and upper limit of the power of 10 the data must breach to convert to scientific notation
            default=mpl.rcParams['axes.formatter.limits']
        preferredOrder: Manually set the power of 10 the data should be converted to. If passed, the limits argument is ignored.

    Returns
        The tuple contains the list of values converted to scientific (in terms of powers of 10) notation
        and the string label to be added to the axis label"""
    # Setting the value of limits if not passed:
    if limits is None:
        limits = mpl.rcParams['axes.formatter.limits']
    label ="" # default value for the axis label
    array = np.array(array) 

    convert = False # Flag to indicate whether data should be converted to scientific notation
    power = 0

    if preferredOrder is None: # Calculating the order based on limits
        # Getting the limits beyond which this algorithm should be applied
        maximum = np.abs(array).max()
        lowerThanRange = maximum < math.sqrt(10) * 10 ** limits[0]
        higherThanRange = maximum > 10 ** limits[1]
        # Debugging:
        # print(f"Sci Notation limits: {limits}\nHigher than range: {higherThanRange}\nLower than range: {lowerThanRange}")
        if lowerThanRange or higherThanRange: # the maximum is out of range
            convert = True
            if maximum == 0:
                print("The passed range cannot be converted to scientific notation since it is only zeroes!")
                convert = False
            power = math.floor(math.log10(maximum))
    else:
        convert = True
        power = preferredOrder

    if convert:
        array = array / (10 ** power) # normalizing with respect to power of 10
        if power != 0:
            label = r"$\times 10^{%i}$ " % (power)
    
    return array, label # unchanged for when the maximum is in range

def ConvertMultipleArraysToScientificNotation(arrayList: List[ArrayLike], limits: None|tuple[int, int] = None, preferredOrder: None|int = None) -> Tuple[List[ArrayLike], str]:
    """Accepts a list of arrays and converts all of them to scientific notation using a common scaling factor (power of 10).

    Args:
        arrayList: List of arrays containing the numerical data to be converted
        limits: Pair (tuple) of integers to set lower and upper limit of the power of 10 the data must breach to convert to scientific notation
            default=mpl.rcParams['axes.formatter.limits']
        preferredOrder: Manually set the power of 10 the data should be converted to. If passed, the limits argument is ignored.

    Returns
        A tuple containing the corresponding list of converted arrays and the string label indicating the scaling factor"""
    arrayLengths = []
    for array in arrayList:
        arrayLengths.append(len(array))
    single_dim_array = np.concatenate(arrayList) # concatenates all arrays into a single array 
    single_dim_array, axisLabel = ConvertToScientificNotation(single_dim_array, limits, preferredOrder)
    # Splitting the array according to the lengths of the original segregation time arrays; returns a list of subarrays as views
    splitIndices = np.cumsum(arrayLengths)[:-1] # array of indices dictating where splitting occurs; cumulative sum of lengths excluding the last sum
    arrayList = np.split(single_dim_array, splitIndices) # just like the original, but converted to scientific notation; list of views
    return arrayList, axisLabel

def ExtractSubstringInParenthesis(string: str) -> str:
    """Extracts the first instance of a substring enclosed in parenthesis. Returns None if none are found."""
    substrings = re.findall("\(.*?\)", string)
    if len(substrings) == 0:
        return None
    else:
        substring = substrings[0][1:-1] # shaving off the parenthesis in the first instance
        return substring
    
def ReadDiameter(numberOfMonomers: int, architecture: str, initializationProcedeure: str) -> float:
    """Reads and returns the diameter corresponding to the architecture passed from the diamaters database"""
    filePath = sysPaths.GetDiameterDatabaseFilePath(numberOfMonomers, initializationProcedeure)
    df = pd.read_csv(filePath)
    architectures = df.iloc[:, 0]
    radiiOfGyration = df.iloc[:, 1]
    diameters = df.iloc[:, 2]
    architectureFound = False
    diameter = 0
    for i in range(len(architectures)):
        if architectures[i] == architecture:
            architectureFound = True
            diameter = diameters[i]
            break
    if not architectureFound:
        print(f"ERROR: The architecture {architecture} was not found in the diameters database!")
        sys.exit(1)
    return diameter


def ReadEffectiveAxisLength(numberOfMonomers: int, architecture: str, cylinderDiameter: float) -> float:
    """
    Reads the effective axis length for a polymer in an infinite cylinder with the given diameter.
    This effeective length is useful to define a length scale on which the segregation criterion can be based.
    Args:
        -numberOfMonomers: The number of monomers present in a single polymer
        -architecture: The name of the topology or architecture
        -cylinderDiameter: The diameter of the cylinder for which the effective length is desired (format: float with two decimal places)
    """
    filePath = sysPaths.GetAxisLengthDatabaseFilePath(numberOfMonomers)
    
    # Reading database file:
    df = pd.read_csv(filePath)
    architectures = df.iloc[: , 0]
    diameters = df.iloc[: , 1]
    halfAxisLengths = df.iloc[: , 2]

    # Searching for combination of architecture and diameter:
    combinationFound = False
    for i in range(len(architectures)):
        if architectures[i] == architecture and math.isclose(diameters[i], cylinderDiameter):
            combinationFound = True
            return 2 * halfAxisLengths[i]
    # If not found:    
    print(f"ERROR: The architecture {architecture} and diameter {cylinderDiameter} combination was not found in the axis lengths database!")
    sys.exit(1)
    return

def GetAxisLength(numberOfMonomers: int, architecture: str, initializationProcedure: str, cylinderDiameter: float|None = None) -> float:
    """
    Returns the axis length for finite cylinders and the effective axis length for infinite cylinders.
    Args:
        -numberOfMonomers: The number of monomers present in a single polymer
        -architecture: The name of the topology or architecture
        -initializationProcedure: The name of the initialization procedure used to generate the initial state
        -cylinderDiameter: (Optional for finite cylinder) The diameter of the cylinder for which the effective length is desired (format: float with two decimal places)
    """
    if sysPaths.IsCylinderInfinite(initializationProcedure):
        if cylinderDiameter is None:
            raise ValueError("cylinderDiameter argument was not passed. cylinderDiameter argument must be passed if the cylinder is infinite.")
        return ReadEffectiveAxisLength(numberOfMonomers, architecture, cylinderDiameter)
    else:
        return segParam.ASPECT_RATIO * ReadDiameter(numberOfMonomers, architecture)


def RescaleArraysByLength(arrayList: list[ArrayLike], length: float) -> list[ArrayLike]:
    """Returns a rescaled version of each array passed in the list.
    The rescaling is done by dividing by the length passed."""
    for i in range(len(arrayList)):
        arrayList[i] = arrayList[i] / length
    return arrayList


def IncludePlotIdentificationLabel(fig: mpl.figure.Figure, label: str, fontsize: int|None = None) -> None:
    """
    Includes a small identification label at the top left corner of the passed figure. 
    This text label is meant to be a label for the plot part of a composite figure in a manuscript.
    Args:
        fig: The matplotlib Figure object on which the label is to be drawn.
        label: The text to be drawn on the figure.
        fontsize: Optional argument for the fontsize of the label text. If unspecified, the value of ID_LABEL_FONTSIZE in Segregation_Parameters.py is taken.
    """

    # loading default fontsize if necessary:
    if fontsize is None:
        fontsize = segParam.ID_LABEL_FONTSIZE

    fig.text(segParam.ID_LABEL_POS[0], segParam.ID_LABEL_POS[1], label, fontsize = fontsize)

def ReadLastLines(filePath: str, N: int = 1) -> list[str]:
    """Reads the last N lines of the file passed and returns them as a tuple"""
    lastLines = []
    # Opening the file in binary mode and seeking to the last line of the file
    with open(filePath, "rb") as file:
        try:
            # points to the last line of the file as long as it contains more than a single line
            file.seek(-2, os.SEEK_END) # goes to two bytes before the end of file
            for i in range(N):
                while file.read(1) != b'\n': # searching if the one byte at the current location is a binary return character
                    file.seek(-2, os.SEEK_CUR) # moving two bytes back from the current position
                # Found a return statement before the desired line:
                lastLines.append(file.readline().decode())
                # Debugging:
                # print(f"ReadLastLines(): Moving {len(lastLines[-1])-10} bytes back.")
                file.seek(-len(lastLines[-1])-10, os.SEEK_CUR) # moving to the beginning of the line stored
        except OSError:
            print(f"Something went wrong while reading the last {N} line(s) of the file {filePath}.")
            print(f"Please make sure the file has atleast {N} lines.")
        lastLines.reverse()
        return lastLines
    
def GetIndexedAxis(axArray: ArrayLike, orderIndex: int, layout: tuple[int, int]) -> mpl.axes.Axes:
    """
    Returns the axis object for the given order index in the specified layout.
    Args:
        axArray (list[mpl.axes.Axes]): A 2D array of axes objects created using the `plt.subplots()` function.
        orderIndex (int): The index of the axis in the flattened array.
        layout (list[list[int]]): The layout of the subplots, where the first element is the number of rows and the second element is the number of columns.

    """
    # Plotting axis indices:
    i = orderIndex // layout[1]
    j = orderIndex % layout[1]
    if layout[0] > 1 and layout [1] > 1: # Both the number of rows and columns are greater than 1
        ax_ij = axArray[i, j]
    else: # Either the number of rows or columns is 1
        ax_ij = axArray[orderIndex]
    return ax_ij
