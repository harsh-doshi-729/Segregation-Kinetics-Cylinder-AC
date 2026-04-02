// This is the config file to from which the copyMixedState.c and createInitialState.c scripts read the values of certain parameters

// copyMixedState.c:

#define USE_CONFIG_PATHS false
#define SOURCE_PATH "/scratch/Harsh/New_Segregation/b200/Previous_Attempts/fbsr_parallel/Arc4/run48/final_configuration.txt"
#define DEST_PATH "/scratch/Harsh/New_Segregation/Check_Concatenation/fbsr_b200/Arc4/post_run48/initial_configuration.txt"
// Some flags:
// Some flags:
#define WRITE_INDIVIDUAL_POLYMER_STATES true
// For most case, both polymers need to copied: set the flag to false
// In the case of copying a state for an R_g / extent run, set the flag to true so that it separates the polymers

#define DO_GHOST_ATOMS_EXIST false
#define DO_INTERPOLYMER_BONDS_EXIST false
#define KEEP_TER_CROSSLINK false

#define ATOM_FLAGS_EXIST true
#define USE_DIFF_BOND_TYPE false // for cross links
#define USE_CROSSLINK_REGIONS false // The cross-linked monomers are used as bounds to define different regions;
// These bounds are added to the existing bounds of 1 and numberOfMonomers

// Sections:
#define VELOCITY_SECTION_EXISTS true
#define MASS_SECTION_EXISTS true

// createInitialState.c:
#define USE_COMMON_TEMPLATE true
// A flag indicating whether the template file to create the
// initial configuration file should be taken from the common b200/ location or from the special_simulation location
