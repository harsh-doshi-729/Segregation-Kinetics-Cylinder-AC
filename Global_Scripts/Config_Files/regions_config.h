// This is a header file that defines the regions used in the region.h header file
#ifndef REGIONS_CONFIG_H // include guard
    #define REGIONS_CONFIG_H

    // defining macros  for the flag variables:
    #define DO_NOT_USE_REGIONS true
    #define DIVIDE_ONLY_FIRST_POLYMER false

    #define ARCHITECTURE "Arc1_10"
    // The region bounds defined as a macro:
    #define REGIONS_ARRAY_LENGTH 2
    #define REGION_BOUNDS_ARRAY {{51, 150}, {151, 50}}
    // #define REGION_BOUNDS_ARRAY {{51, 150}, {1, 10}, {11, 20}, {21, 30}, {31, 40}, {41, 50}, {151, 160}, {161, 170}, {171, 180}, {181, 190}, {191, 200}}
    // #define REGION_BOUNDS_ARRAY {{51, 62}, {88, 112}, {138, 150}, {1, 50}, {63, 87}, {113, 137}, {150, 200}}
    // This is the array for the region bounds within a polymer. Both bounds are included in the region
    // Note: The order of the indices in the pair of bounds matters. The 1st of the pair is the lower index and the second is the upper index. 

    #define NON_LOOP_REGIONS 3 // The number of non-looped regions in the above array; all of them occur at the start of the array
#endif