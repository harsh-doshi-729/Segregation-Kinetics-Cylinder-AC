# This is a python script that will be used as a config file to import the values of the parameters 
# used in scripts
import system_file_paths as sysPaths
import cluster_file_paths as clusPaths

NUMBER_OF_RUNS = 50
# Aspect ratio of cylindrical confinement:
ASPECT_RATIO = 5 # Ratio of cylinder length to diameter

# Parameters for segregation time criterion:
# fraction of the box length:
FIRST_THRESHOLD = 0.45
SECOND_THRESHOLD = 0.40
THIRD_THRESHOLD = 0 # ignoring remixing effects
INTERVAL_LENGTH = 0.50 # The fraction of the total run over which the average of CoM distance should be calculated to compare to second threshold
CRITERION_STRING = "f045_s040_t00" # this must match the above threshold values and interval length;
# If L = 0.5; then omit the '_Lxx' from the string

# Parameters for segregation in an infinite cylinder:
INF_F_THRESHOLD = 0
INF_S_THRESHOLD = 0.01
INF_CRITERION_STRING = "inf_f000_s001"

# CoM Time series plotting style:
SPLIT_COM_DISTANCE = False # If true, the polymer CoM time series is plot for each polymer separately in a single plot
# If false, the CoM distance is plotted as a function of time

# Tau_0 relaxation time of one monomer:
TAU_0 = 100 # approximate relaxation time in terms of iterations of Langevin Dynamics with timestep 0.01

# Architectures that count as or are based on linear polymers:
LINEAR_ARCHITECTURES = ["Arc_Lin"]

# Architectures of interest: to perform operations on a few select architectures
AOI_LIST = ["Arc0", "ArcI-8", "Arc2", "Arc_Loop_5", "Arc1_10"]
# AOI_LIST = ["Arc0", "Arc_Loop_5"]
# AOI_LIST = ["Arc0", "Arc1_1", "ArcR8_1_9", "Arc4", "ArcI-8"]
# AOI_LIST = ["Arc0", "ArcR8_1_19", "ArcR8_1_9", "ArcR8_3_7", "ArcI-8"]
# AOI_LIST = ["Arc0", "Arc2", "Arc1_2_2", "ArcCW_6_7_7"]
# AOI_LIST = ["Arc0", "Arc1_1", "ArcML_2_20", "ArcML_5_20", "Arc1_10"]
# AOI_LIST = ["Arc0", "Arc1_1", "ArcR8_1_9", "Arc4", "ArcI-8", "Arc2", "Arc_Loop_5", "Arc1_10"]
# AOI_LIST = ["6.24_Arc0", "6.24_ArcI-8", "6.24_Arc2", "6.24_Arc_Loop_5", "6.24_Arc1_10"]
LOOP_MODIFIER = "multi_loops" # "single_loop" # The file name modifier while saving box plots of a particular set of architectures

SORT_ARCHITECTURES = True # Flag to indicate whether the AOI_LIST should be sorted according to their names

# Architecture aliases: A dictionary to store the actual names and the alias names of the architectures
# The keys: the names given to give the runs and label directories
# The values are the (hopefully) polished/sophisticated names or aliases that are to be reported
# ALIASES = {"Arc1_1": "ArcR8-1-19", "ArcR8_1_19": "ArcR8-1-19", "ArcR8_1_9": "ArcR8-1-9", "Arc4": "ArcR8-3-7", "ArcR8_3_7": "ArcR8-3-7", "ArcI-8": "ArcR8-1-1"} # Old nomenclature
# ALIASES = {"Arc1_1": "Arc-1-1", "ArcR8_1_19": "Arc-1-1", "ArcR8_1_9": "Arc-1-1", "Arc4": "Arc-1-1", "ArcR8_3_7": "Arc-1-1", "ArcI-8": "Arc-1-1"}
ALIASES = {"Arc0": "Arc-0", "ArcI-8": "Arc-1-1", "Arc2": "Arc-1-2", "Arc_Loop_5": "Arc-1-5", "Arc1_10": "Arc-1-10"}
# The internal loop size breakdown for all architectures:
ARC_BREAKDOWN_200 = {"Arc0": "[200]", "Arc1_1": "[190-10]", "ArcR8_1_19": "[190-10]", "ArcR8_1_9": "[180-20]",
                     "Arc4": "[140-60]", "ArcR8_3_7": "[140-60]", "ArcI-8": "[100-100]",
                     "Arc2": "[100-50]", "Arc_Loop_5": "[100-20]", "Arc1_10": "[100-10]"}

ARC_BREAKDOWN_500 = {"Arc0": "[500]", "Arc1_1": "[475-25]", "ArcR8_1_19": "[475-25]", "ArcR8_1_9": "[450-50]",
                     "Arc4": "[350-150]", "ArcR8_3_7": "[350-150]", "ArcI-8": "[250-250]",
                     "Arc2": "[250-125]", "Arc_Loop_5": "[250-50]", "Arc1_10": "[250-25]"}

def GetAliasArchitecture(architecture: str) -> str:
    """Returns the corresponding alias architecture to the passed architecture if it exists in ALIASES"""
    if architecture in ALIASES.keys():
        return ALIASES[architecture]
    else:
        return architecture
    
def GetAliasArchitectureList(architectures: list[str]) -> list[str]:
    """Returns a list of architectures containing the corresponding alias architectures to the ones in the passed list"""
    aliasArchitectures = []
    for architecture in architectures:
        aliasArchitectures.append(GetAliasArchitecture(architecture))
    return aliasArchitectures

def GetPaperFiguresTitle(numberOfMonomers: int, architecture: str, special_simulation: str = sysPaths.SPECIAL_SIMULATION) -> str:
    """Returns the general minimalistic title for paper-ready plots. 
    The title incorporates the number of monomers, architecture name, and special simulation."""
    breakdown = ""
    if numberOfMonomers == 200:
        breakdown = f" {ARC_BREAKDOWN_200[architecture]}"
    elif numberOfMonomers == 500:
        breakdown = f" {ARC_BREAKDOWN_500[architecture]}"
    return f"{GetAliasArchitecture(architecture)}{breakdown} - {GetAliasSimulation(special_simulation)}"
# TODO: Add number of runs


# The special simulations used in comparing box plots:
# SPECIAL_SIMULATIONS = ["fene_recenter", "fene_LJ", "fene_fbsr", "replication_like"]
SPECIAL_SIMULATIONS = ["fene_fbsr", "fbsr_antiparallel"]
# AOI_COMPARE = ["Arc0", "ArcI-8", "Arc2", "Arc_Loop_5", "Arc1_10"] # An arc list for each special simulation
PLOT_COMPARISON = False # A flag to indicate whether to plot the comparison for different architectures/runs
AOI_COMPARE = ["Arc0", "ArcI-8", "Arc2", "Arc_Loop_5", "Arc1_10"] # ["Arc1_10", "Arc_Loop_5", "Arc2", "ArcI-8", "Arc0"]
# Custom folder paths for any comparison scripts: as a dictionary of <identifier>-<custom path> key-value pair
COMPARE_CUSTOM_PATHS = {"Arc-1-5": "/home/harsh/ComputationalPhysics/Polymer_Physics/LAMMPS_runs/cluster-data/new_segregation/b200/Previous_Attempts/inf_recenter/Arc_Loop_5/run7/",
                        "Arc-1-10": "/home/harsh/ComputationalPhysics/Polymer_Physics/LAMMPS_runs/cluster-data/new_segregation/b200/Previous_Attempts/inf_recenter/Arc1_10/run1/"}

ALIAS_SPECIAL_SIMULATIONS = {"fene_recenter": "Recenter-COM", "fene_LJ": "Mutual-Attraction", "fene_fbsr": "Glued-Monomers", "replication_like": "Ladder-like", "inf_recenter": "Recenter-COM", "inf_LJ": "Mutual-Attraction", "inf_fbsr": "Glued-Monomers", "inf_replication": "Ladder-like"}

def GetAliasSimulation(special_simulation: str = sysPaths.SPECIAL_SIMULATION) -> str:
    """Returns the alias of the special simulation passed.
    The default value of the special simulation is the one set in the system_file_paths.py file.
    If the special simulation does not have a set alias, the passed name itself is returned."""
    if special_simulation in ALIAS_SPECIAL_SIMULATIONS.keys():
        return ALIAS_SPECIAL_SIMULATIONS[special_simulation]
    else:
        return special_simulation
# The number of rows and columns for the comparison plot of segregation times for different architectures:
# Or the number of rows and columns for the comparison box plots of different special simulations
NROWS = 1
NCOLS = 3

# Groups of different runs:
USE_GROUPS = False
NGROUPS = 2 # number of groups
# List of lists of indices; sublists contains run indices belonging to the same group
# GROUP_INDICES = [[1, 2, 3, 8, 9, 13, 14, 15, 16, 17, 20, 22, 23, 24, 28, 30, 32, 33, 35, 36, 38, 39, 40, 41, 42, 45, 48]] # Arc2
# GROUP_INDICES = [[2, 3, 4, 7, 9, 10, 11, 13, 14, 17, 19, 21, 22, 23, 24, 25, 27, 28, 38, 40, 42, 43, 44, 45, 49, 50]] # Arc2, fene_recenter (200)
# GROUP_INDICES = [[4, 6, 9, 10, 12, 16, 19, 20, 23, 24, 25, 28, 29, 35, 39, 44, 46]] # Arc_Loop_5, fene_recenter (200)
# GROUP_INDICES = [[1, 2, 3, 6, 12, 16, 19, 20, 25, 27, 28, 31, 33, 35, 36, 37, 39, 40, 41, 45, 46, 47, 50]] # Arc2, fene_recenter (500)
GROUP_INDICES = [[4, 11, 13, 15, 17, 19, 22, 26, 27, 31, 33, 37, 45, 47]] # Arc_Loop_5, fene_recenter (500)
# Contains NGROUP-1 lists, the last list is the implicitly rest of the indices
GROUP_LABELS = ["Parallel", "Antiparallel"] # Corresponding labels for the groups

def UpdateGroupIndices(universe: set) -> None:
    """Updates the defined GROUP_INDICES to add the complement indices from the passed universe set"""
    if len(GROUP_INDICES) < NGROUPS:
        indexSet = set()
        for indices in GROUP_INDICES:
            indexSet = indexSet | set(indices) # union of sets
        GROUP_INDICES.append(list(universe - indexSet))

# Using custom labels for the Architecture segregation times box/mean plots
USE_CUSTOM_LABELS = True
CUSTOM_AXIS_LABEL = "Number of small loops"
# CUSTOM_AXIS_LABEL = "Number of monomers in small loop"
CUSTOM_LABELS = {"Arc0": 0, "ArcI-8": 1, "Arc1_1": 1, "ArcR8_1_19": 1, "Arc2": 2, "ArcML_2_20": 2, "Arc_Loop_5": 5, "ArcML_5_20": 5, "Arc1_10": 10}
# CUSTOM_LABELS = {"Arc0": 0, "Arc1_1": 10, "ArcR8_1_19": 10, "ArcR8_1_9": 20, "Arc4": 60, "ArcR8_3_7": 60, "ArcI-8": 100}
# CUSTOM_LABELS = {"Arc0": 0, "Arc1_1": 25, "ArcR8_1_19": 25, "ArcR8_1_9": 50, "Arc4": 150, "ArcR8_3_7": 150, "ArcI-8": 250}

def GetCustomLabel(architecture: str) -> str:
    """Returns the custom label corresponding to the passed architecture"""
    if architecture in CUSTOM_LABELS.keys():
        return CUSTOM_LABELS[architecture]
    else:
        raise Exception(f"The architecture {architecture} does not exist in the custom labels dictionary!")
    
def GetCustomLabels(architectures: list[str]) -> str:
    """Returns a list of labels corresponding to the list of architectures passed"""
    customLabels = []
    for architecture in architectures:
        customLabels.append(GetCustomLabel(architecture))
    return customLabels

# Flag to include flier information in box plots:
INCLUDE_FLIERS_INFORMATION = True

# Plot inset of segregation trajectory comparisons:
PLOT_INSET = True
INSET_POSITION = [0.4, 0.1, 0.5, 0.5] # [x0, y0, width, height]

# Plot only Trajectory comparison using AOI_COMPARE:
ONLY_TRAJECTORY_COMPARISON = True

# Whether to use loglog plots for MSD:
USE_LOGLOG_PLOT = False # Takes precedence over semilog flag
USE_SEMILOG_PLOT = True

# Whether all the individual runs should be plotted in the background as light grey lines (while calculating the mean squared CoM displacement):
PLOT_ALL_RUNS = True

# The mode of the Squared COM trajectory: Distance or Displacement
# Distance: (z(t))^2
# Displacement: (z(t) - z(0))^2
COM_MODE = "Distance" # "Distance" or "Displacement"

# Whether to normalize the Std. Dev. of Squared COM by the mean squared COM displacement
NORMALIZE_STD_DEV = False

# Whether to calculate or read the COM distance distribution from file:
READ_COM_DISTANCE_DISTRIBUTION = False
# Whether to thw write the calculate COM distance distribution to file:
WRITE_COM_DISTANCE_DISTRIBUTION = True

# A time (in no of iterations) cutoff introduced to excluded large contributions from the segregated state in the COM distance distribution
USE_DISTRIBUTION_CUTOFF = True
COM_DISTRIBUTION_CUTOFF = 2 * 10 ** 7 # 3 * 10 ** 6; for 200

# A set of common x limits and y limits for any plot:
USE_COMMON_AXIS_LIMITS = True
X_LIMS = [-0.55, 0.55]
Y_LIMS = [-0.05, 1.5]

# Preferred order to convert to scientific notation:
PREFERRED_ORDER = 4

# Free Energy parameters: A dictionary such that the key is the file label and the value is the label that is to be displayed
# FREE_ENERGY_LABELS = {"region_CoM_1_3": "Loop 1 and 3", "region_CoM_2_4": "Loop 2 and 4"}
FREE_ENERGY_LABELS = {"region_CoM_1_3": "Small Loops", "region_CoM_2_4": "Big Loops"}
# FREE_ENERGY_LABELS = {"region_CoM_1_2": "z CoM Distance"}
# A dictionary such that the key is the free energy label prefix, and the value is the axis label for the relevant reaction coordinate
FREE_ENERGY_AXIS_LABELS = {"region_CoM": r"$\Delta z_{COM}$", "region_overlap": r"Overlap distance $z$"}

# Parameter for plotting monomer distributions for initial mixed states:
USE_SINGLE_SNAPSHOT = False # Set this flag to true if the regional monomer distributions from a single snapshot need to be plotted
# The above distribution is plotted after reading from a file where the single snapshot distribution is already computed
# This feature now has been superceded by plotSingleSnapshotDistribution.py. The flag should be False.
# Note: The old computed distribution may not match with the distribution computed by new script since the
#       input positions might be different for the old method
PLOT_SINGLE_SNAPSHOT_DENSITY = False # A flag to indicate whether to plot the probability distribution (if True) or a histogram (if False)

# For plotting graphs with rescaled lengths : normalized with respect to the box length
if sysPaths.SPECIAL_SIMULATION[0:4] == "inf_":
    RESCALE_LENGTHS = False # Rescale distances by the cylinder length
else:
    RESCALE_LENGTHS = True
# Set false for infinite cylinder
RESCALE_RADIUS = True # Rescale radial coordinates by the cylinder radius

# File extension for saving figures and plots:
FIG_EXT = ".png" # ".eps"
SHOW_TITLE = False # The title might not be relevant for including plots in a paper or poster
USE_MARGINS = False
PLOT_MARGINS = {'left': 0.14, 'bottom': 0.14, 'right': 0.96, 'top': 0.96} # What for?
# PLOT_MARGINS = {'left': 0.03, 'bottom': 0.14, 'right': 0.96, 'top': 0.96} # For monomer density

# Custom commands for writing an SFTP script; currently for backing up New_Segregation zip files
SFTP_S_BASE_FOLDER = "/scratch/Harsh/New_Segregation/" # Base NEW_SEGREGATION folder on the 175 cluster
SFTP_LOCAL_BASE_FOLDER = "/media/soft-matter-group/Expansion/Harsh/new_segregation/" # External drive
def Get_SFTP_Commands(runIndex: int) -> str:
    """Returns the set of SFTP commands to be used for a single run"""
    return f"""
    lmkdir run{runIndex}/
    lcd run{runIndex}/
    cd run{runIndex}/
    get visual.dump.zip
    lcd ..
    cd ..
    """

# Special simulations and architectures to be backed up locally:
BACKUP_SPECIAL_SIMULATIONS = ["fene_fbsr", "replication_like"]
BACKUP_ARCHITECTURES = ["Arc0", "ArcI-8", "Arc2", "Arc_Loop_5", "Arc1_10"]

# Parameters to include small label in figure:
INCLUDE_PLOT_IDENTIFICATION_LABEL = True
ID_LABEL_POS = (0.005, 0.93) # For subplot figures: (0.01, 0.95) # In Figure units from (0,0) to (1,1)
ID_LABEL_TEXT = "(c)"
ID_LABEL_FONTSIZE = 36 # for subplot figures: 20
