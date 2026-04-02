// This is the config file to from which the copyMixedState.c script reads the value of certain parameters

#define USE_CONFIG_PATHS true
#define SOURCE_PATH "/scratch/Harsh/New_Segregation/b200/Previous_Attempts/fbsr_parallel/Arc4/run48/final_configuration.txt"
#define DEST_PATH "/scratch/Harsh/New_Segregation/Check_Concatenation/fbsr_b200/Arc4/post_run48/initial_configuration.txt"
// Some flags:
#define DO_GHOST_ATOMS_EXIST false
#define DO_INTERPOLYMER_BONDS_EXIST false
#define KEEP_TER_CROSSLINK false

#define ATOM_FLAGS_EXIST true
#define USE_DIFF_BOND_TYPE false // for cross links
#define USE_CROSSLINK_REGIONS true
// Sections:
#define VELOCITY_SECTION_EXISTS true
#define MASS_SECTION_EXISTS true