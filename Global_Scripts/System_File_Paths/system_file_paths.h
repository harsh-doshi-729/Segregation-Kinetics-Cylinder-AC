// This header files defines all the different local filepaths as macros for my scripts to use.
// This will make switching between different systems seamless

#ifndef SYSTEM_FILE_PATHS_H
    #define SYSTEM_FILE_PATHS_H // include guard

    // IMPORTANT NOTE: Please set the value of the macro below to the path of the directory where the project is located on your system. This will make it so that all the file paths defined in this header file are correct for your system. For example, if the project is located in "/home/harsh/Segregation-Kinetics-Cylinder-AC/", then set BASE_DIR to "/home/harsh/Segregation-Kinetics-Cylinder-AC/".
    #define BASE_DIR "/home/harsh/Segregation-Kinetics-Cylinder-AC/" // The base directory where the project is located; all other file paths are defined based
    

    // Utility files:
    // #define RANDOM "/home/harsh/ComputationalPhysics/Week_1_Assignment/random.h"
    // NOTE: Please ensure to add the option -I <path/to/project/directory> while compiling any of the C scripts that import this header file
    // This makes it so that the compiler looks for the header files in the project directory as well, and thus, can find the header file paths defined in this file
    #define GENERAL_FILE_METHODS "/Global_Scripts/Utility/GeneralFileMethods.h"
    #define LAMMPS_POSITION_FILE "Global_Scripts/Utility/LAMMPS_Position_File.h"
    #define MONOMER_DISTRIBUTION "Global_Scripts/Utility/MonomerDistribution.h"
    #define CROSS_LINKS_DATABASE "Global_Scripts/Utility/Cross-Links_Database.h"
    #define LOOP_CUTTING "Global_Scripts/Utility/LoopCutting/Scripts/loop_cutting.h"
    #define SECTION "Global_Scripts/Utility/Section.h"
    #define MODULO "Global_Scripts/Utility/modulo.h"
    #define REGION "Global_Scripts/Utility/region.h"
    #define SEGREGATION_TIMES "Global_Scripts/Utility/Segregation_Times_File.h"
    #define EXTENT "Global_Scripts/Utility/extent.h"
    #define PAIR_CORRELATION "Global_Scripts/Utility/pair_correlation.h"

    // Config files:
    #define REGION_CONFIG "Global_Scripts/Config_Files/regions_config.h"
    #define COPY_CONFIG "Global_Scripts/Config_Files/copy_config.h"
    #define STATE_CONFIG "Global_Scripts/Config_Files/state_config.h"
    
    #define SPECIAL_SIMULATION "fene_recenter" // Optional macro for special simulations; Set it to blank if no need for special simulation folder
    #define R_G_SIMULATION "inf_recenter"
    
    // Important directories as strings:
    #define HOME_DIR "/home/harsh/Segregation-Kinetics-Cylinder-AC/"
    #define SEGREGATION "Segregation/"
    #define CHECK_CONCATENATION "Check_Concatenation/"
    #define CREATE_INITIAL_STATES "Create_Initial_States/"



    #define ANALYSIS "Global_Scripts/Analysis/" // This Analysis folder can be any of the Analysis folder in Linear_Chains on the cluster

    // Segregation criterion:
    #define SEGREGATION_CRITERION "f045_s040_t00"

#endif
