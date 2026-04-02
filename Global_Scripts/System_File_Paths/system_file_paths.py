# This Python script defines constant variables for the local paths that other Python script may use

SPECIAL_SIMULATION="inf_recenter" # The folder in which the simulation data files lie; it is not blank when a special simulation's data needs to be accessed
R_G_SIMULATION="" # The folder in which the simulation data files lie; it is not blank when a special simulation's data needs to be accessed

NEW_SEGREGATION = f"/home/harsh/ComputationalPhysics/Polymer_Physics/LAMMPS_runs/cluster-data/new_segregation/"
CREATE_INITIAL_STATES = f"/home/harsh/ComputationalPhysics/Polymer_Physics/LAMMPS_runs/cluster-data/new_segregation/Create_Initial_States/"
CHECK_CONCATENATION = f"/home/harsh/ComputationalPhysics/Polymer_Physics/LAMMPS_runs/cluster-data/new_segregation/Check_Concatenation/"
R_G = "/home/harsh/ComputationalPhysics/Polymer_Physics/LAMMPS_runs/cluster-data/new_segregation/R_g/"

LOOP_CUT_INITIAL_STATES = "/home/harsh/ComputationalPhysics/Polymer_Physics/LAMMPS_runs/cluster-data/LoopCutting/Initial_States/"
CUT_STATES = "/home/harsh/ComputationalPhysics/Polymer_Physics/LAMMPS_runs/cluster-data/LoopCutting/Cut_States/"

CLUSTER_INITIALIZATION = "/home/harsh/ComputationalPhysics/Polymer_Physics/LAMMPS_runs/cluster-data/initialization/"
INITIALIZATION = "/home/harsh/ComputationalPhysics/Polymer_Physics/LAMMPS_runs/initialization/"

OUTPUT_FILES = "/home/harsh/ComputationalPhysics/Polymer_Physics/Output_Files/"
LAMMPS_RUNS = "/home/harsh/ComputationalPhysics/Polymer_Physics/LAMMPS_runs/"

# The path to the folder where the git repo is based
POLYMER_PHYSICS = "/home/harsh/ComputationalPhysics/Polymer_Physics/"

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
        return f"{baseFolder}b{numberOfMonomers}/Previous_Attempts/{simulation_name}/"

def GetDiameterDatabaseFilePath(numberOfMonomers: int) -> str:
    return f"{NEW_SEGREGATION}b{numberOfMonomers}/Diameters.csv"

def IsCylinderInfinite(special_simulation: str|None = None) -> bool:
    """Returns whether the passed special_simulation corresponds to the case of an infinite cylinder.
    If no argument is passed, the default SPECIAL_SIMULATION value is assumed."""
    if special_simulation == None:
        special_simulation = SPECIAL_SIMULATION
    # If the special simulation has the 'inf_' prefix, then it corresponds to an infinite cylinder case
    return (special_simulation[0: 4] == "inf_")

def GetAxisLengthDatabaseFilePath(numberOfMonomers: int) -> str:
    return f"{NEW_SEGREGATION}b{numberOfMonomers}/AxisLengths.csv"
