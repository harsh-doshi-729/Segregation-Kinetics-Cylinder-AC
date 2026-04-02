# This is a config file that will be useful for plotting the monomer densities  or CoM timesries for polymers with multiple regions

import sys
from system_file_paths import POLYMER_PHYSICS


USE_REGIONS = True # A flag to indicate whether the polymers are subdivided into regionss
NUMBER_OF_POLYMER_REGIONS = 2 # The number of regions in a single polymer
NUMBER_OF_REGIONS = 4 # The total number of distinct regions in the entire system
if USE_REGIONS:
    REGION_LABELS = ["Large Loop - P1", "Small Loops - P1", "Large Loop - P2", "Small Loops - P2"]
    # REGION_LABELS = ["Large Loop - P1", "L1 - P1", "L2 - P1", "L3 - P1", "L4 - P1", "L5 - P1",
                    #   "Large Loop - P2", "L1 - P2", "L2 - P2", "L3 - P2", "L4 - P2", "L5 - P2"]
    # REGION_LABELS = ["Large Loop - P1", "Small Loop - P1", "Large Loop - P2", "Small Loop - P2"]
    # REGION_LABELS = ["Large Loop - P1", "Small Loop 1 - P1", "Small Loop 2 - P1", "Large Loop - P2", "Small Loop 1 - P2", "Small Loop 2 - P2"] # Arc2 verbose
    # REGION_LABELS = ["Large Loop - P1", "Small Loop (1-50) - P1", "Small Loop (150-200) - P1", "Large Loop - P2", "Small Loop (1-50) - P2", "Small Loop (150-200) - P2"] # Arc2 like; for other arcs
    # REGION_LABELS = ["Loop 1 - P1", "Loop 2 - P1", "Loop 1 - P2", "Loop 2 - P2"]
    # REGION_LABELS = ["Linker 1 - P1", "Linker 2 - P1", "Small Loops 1 - P1", "Small Loops 2 - P1", "Linker 1 - P2", "Linker 2 - P2", "Small Loops 1 - P2", "Small Loops 2 - P2"]
    # REGION_LABELS = ["Linker 1 - P1", "Linker 2 - P1", "Linker 3 - P1", "Large Loops 1 - P1", "Small Loop 1 - P1", "Small Loop 2 - P1", "Linker 1 - P2", "Linker 2 - P2", "Linker 3 - P2", "Large Loops - P2", "Small Loop 1 - P2", "Small Loop 2 - P2"]
else:
    NUMBER_OF_REGIONS = 2
    NUMBER_OF_POLYMER_REGIONS = 1
    REGION_LABELS = ["Polymer 1", "Polymer 2"]
REGION_BOUNDS = [[51, 150], [151, 50]] # The monomer bounds of the regions in a single polymer
# REGION_BOUNDS = [[126, 375], [376, 125]] # For 500 monomers
# REGION_BOUNDS = [[126, 375], [376, 425], [426, 475], [476, 25], [26, 75], [76, 125]] # For 500 monomers
# REGION_BOUNDS = [[51, 150], [1, 50], [151, 200]]
# REGION_BOUNDS = [[11, 190], [191, 10]]
# REGION_BOUNDS = [[26, 175], [176, 25]]
# REGION_BOUNDS = [[81, 19], [20, 80]]
# REGION_BOUNDS = [[6, 195], [196, 5]]
# REGION_BOUNDS = [[31, 20], [21, 30]]  
# REGION_BOUNDS = [[13, 488], [489, 12]]
# REGION_BOUNDS = [[26, 475], [476, 25]]
# REGION_BOUNDS = [[76, 425], [426, 75]]
# REGION_BOUNDS = [[62, 188], [313, 438], [439, 61], [189, 312]]
# REGION_BOUNDS = [[88, 162], [338, 412], [413, 87], [163, 337]]
# REGION_BOUNDS = [[126, 185], [219, 281], [345, 375], [375, 125], [186, 218], [282, 344]]


def GetRegionID(monomerIndex: int, numberOfMonomers: int) -> int:
    """Returns the region ID according to the region definitions by the regions array; indexing starts from 0
    The monomer index is expected to start from 1"""

    polymerIndex = int((monomerIndex - 1) // numberOfMonomers)
    reducedMonomerID = (monomerIndex - 1) % numberOfMonomers + 1
    if not USE_REGIONS:
        return polymerIndex # Polymer as an entire region
    
    regionIndex = -1
    for i in range(NUMBER_OF_POLYMER_REGIONS):
        if REGION_BOUNDS[i][0] <= REGION_BOUNDS[i][1]: # Normal (continuous) case
            # checking if monomer lies between the bounds:
            if reducedMonomerID >= REGION_BOUNDS[i][0] and reducedMonomerID <= REGION_BOUNDS[i][1]:
                regionIndex = i
                break
        else: # Abnormal case: The discontinuity if monomer index occurs in this region
            # checking if monomer lies beyond the bounds:
            if reducedMonomerID >= REGION_BOUNDS[i][0] or reducedMonomerID <= REGION_BOUNDS[i][1]:
                regionIndex = i
                break
    if regionIndex == -1:
        print(f"ERROR: The region ID for monomer {monomerIndex} could not be found! Please check the region bounds in the config file. Terminating.")
        sys.exit(1)
    else:
        return int(NUMBER_OF_POLYMER_REGIONS * polymerIndex + regionIndex)

