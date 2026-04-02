import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

filePath = ""

def SetPath():
    """Sets the path to the closestDistance.txt file by accepting input from the command line"""
    global filePath
    if len(sys.argv) < 2:
        print("Not enough arguments specified while invoking the script! Please pass the relative path to the closestDistance.txt file in the following format:")
        print("python /<path>/<this_script.py> <relative_path>/closestDistance.txt")
        quit() # terminating the script
    
    filePath = sys.argv[1]

def PlotClosestDistance(timeSteps: list[int], distances: list[float]):
    """Plots and displays the graph of the closest distances vs timesteps"""
    fig, ax = plt.subplots()
    ax.plot(timeSteps, distances, marker = '.')
    ax.set_title("Closest Distance between any two monomers in the simulation")
    ax.set_ylabel("Closest Distance (LJ Units)")
    ax.set_xlabel("Timesteps")

    plt.show(block = True)

def ReadData() -> tuple[list[int], list[float]]:
    """Reads and returns the tuple of two arrays: (timesteps, closestDistances) from the inputted path"""
    df = pd.read_csv(filePath)

    timeSteps = df.iloc[:, 0]
    distances = df.iloc[: ,1]

    return (timeSteps, distances)

def __main__():
    """Main Method"""
    SetPath()
    timeSteps, distances = ReadData()
    PlotClosestDistance(timeSteps, distances)

if __name__ == '__main__':
    __main__()