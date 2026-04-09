# This Python script defines constant variables for the local paths that other Python script may use

SPECIAL_SIMULATION="inf_recenter" # The folder in which the simulation data files lie; it is not blank when a special simulation's data needs to be accessed
R_G_SIMULATION="" # The folder in which the simulation data files lie; it is not blank when a special simulation's data needs to be accessed

import os
import sys
BASE_DIR = os.environ.get("SEG_BASE_DIR") # Reading the base directory from the environment variable;
if BASE_DIR is None or BASE_DIR == "":
    print("SEG_BASE_DIR not set! Please set the environment variable SEG_BASE_DIR to the absolute path of the project root directory.")
    print("For eg.: $ export SEG_BASE_DIR=/path/to/Segregation-Kinetics-Cylinder-AC/")
    sys.exit(1)
SEGREGATION = f"{BASE_DIR}/Segregation/"
CREATE_INITIAL_STATES = f"{BASE_DIR}/Create_Initial_States/"
CHECK_CONCATENATION = f"{BASE_DIR}/Check_Concatenation/"
GLOBAL_SCRIPTS = f"{BASE_DIR}/Global_Scripts/"
R_G = f"{BASE_DIR}/R_g/"

LOOP_CUT_INITIAL_STATES = f"{BASE_DIR}/LoopCutting/Initial_States/"
CUT_STATES = f"{BASE_DIR}/LoopCutting/Cut_States/"

CLUSTER_INITIALIZATION = f"{BASE_DIR}/initialization/"
INITIALIZATION = f"{BASE_DIR}/initialization/"

OUTPUT_FILES = f"{BASE_DIR}/Output_Files/"
LAMMPS_RUNS = f"{BASE_DIR}/LAMMPS_runs/"

# Folder to store backup files (zip files):
BACKUP_NEW_SEGREGATION = "/home/harsh/ComputationalPhysics/Backup_Polymer_Physics/New_Segregation/"
BACKUP_CIS = "/home/harsh/ComputationalPhysics/Backup_Polymer_Physics/Create_Initial_States/"

def GetFolder(baseFolder: str, numberOfMonomers: int, simulation_name: str = "") -> str:
    """Returns the path of the special simulation folder located in Previous_Attempts if simulation_name is passed.
    Otherwise returns the base folder location with the number of monomers
    Args:
        -baseFolder: Path to a folder that contains the system heirarchy of b<N>/Previous_Attempts/... For example: sysPaths.NEW_SEGREGATION
        -numberOfMonomers: N or the number of monomers in a single polymer
        -simulation_name: The name for the special simulation (optional)
    """
    if len(simulation_name) == 0: # No special simulation: default location
        return f"{baseFolder}b{numberOfMonomers}/"
    else:
        return f"{baseFolder}b{numberOfMonomers}/{simulation_name}/"

def GetDiameterDatabaseFilePath(numberOfMonomers: int, initializationProcedure: str) -> str:
    diametersDir = f"{SEGREGATION}b{numberOfMonomers}/"
    if IsCylinderInfinite(initializationProcedure):
        return f"{diametersDir}c_Diameters.csv" # Diameters for constant confinement (infinite cylinder)
    else:
        return f"{diametersDir}Diameters.csv" # Diameter for constant volume fraction (finite cylinder)

def IsCylinderInfinite(special_simulation: str|None = None) -> bool:
    """Returns whether the passed special_simulation corresponds to the case of an infinite cylinder.
    If no argument is passed, the default SPECIAL_SIMULATION value is assumed."""
    if special_simulation == None:
        special_simulation = SPECIAL_SIMULATION
    # If the special simulation has the 'inf_' prefix, then it corresponds to an infinite cylinder case
    return (special_simulation[0: 4] == "inf_")

def GetAxisLengthDatabaseFilePath(numberOfMonomers: int) -> str:
    return f"{SEGREGATION}b{numberOfMonomers}/AxisLengths.csv"
