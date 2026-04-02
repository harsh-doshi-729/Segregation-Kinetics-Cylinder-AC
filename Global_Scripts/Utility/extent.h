// This header file contains variables and functions to aid in calculating the spatial extent of two polymers
// This is useful in discerning whether two polymers are overlapped or segregated

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifndef EXTENT_H // include guard
    #define EXTENT_H
    // Importing the region bounds from a config file:
    #include "../Config_Files/regions_config.h"
    #include "../System_File_Paths/system_file_paths.h"
    #include MODULO
    #include REGION
    #include LAMMPS_POSITION_FILE

    #define OVERLAP_THRESHOLD 1 // Threshold set to sigma; there must be atleast a monomer distance between regions 

    // A struct to hold the minimum and maximum coordinates belonging to a region
    struct Pair {
        double min;
        double max;
    };

    typedef struct Pair Pair;

    // Prints the values of pair to console
    void PrintPair(Pair pair)
    {
        printf("Pair min = %lf\nPair max = %lf\n", pair.min, pair.max);
    }

    // Sets the values of the pair of min and max values for each region for each spatial dimension from the passed simulation_data
    // Arguments: 
    // -2D Array of dim (numberOfRegions * 3) for three spatial directions of each region
    // -the region_data object with the region information
    // -the simulation_data object with the position data of the entire system
    void SetMinMax(Pair pairs[][3], region_data regionData, simulation_data read_data)
    {
        // pairs is a 2D array with numberOfRegions * 3 elements
        // The first index is the region index, and the second is the spatial dimension

        // initializing pairs:
        for(int n = 0; n < regionData.totalNumberOfRegions; n++)
        {
            for(int j = 0; j < 3; j++)
            {
                pairs[n][j].min = __DBL_MAX__;
                pairs[n][j].max = -__DBL_MAX__;
            }
        }

        //iterating over all positions:
        for(int i = 0; i < read_data.numberOfMonomers; i++)
        {
            int regionID = GetRegionID(regionData, i+1);
            // checking the minimum and maximum for that particular region
            for(int j = 0; j < 3; j++)
            {
                double position = read_data.bead_positions[i][j]; // For easier readability
                if(position < pairs[regionID][j].min)
                    pairs[regionID][j].min = position;
                if(position > pairs[regionID][j].max)
                    pairs[regionID][j].max = position;
            }
        }
        // Pairs should be set for all regions!
    }

    // Returns the minimum distance between monomers from Pairs of two regions
    // If the regions do not overlap, the minimum distance is simply the distance between the closest monomers
    // If the regions overlap, then the distance is negative and its magnitude indicates the length of the overlapping distance
    double GetMinimumDistance(Pair pair1, Pair pair2)
    {
        Pair lowerPair, higherPair; // Pairs assigned based on the relative locations of the two regions

        // Assigning location:
        if(pair1.min < pair2.min)
        {
            lowerPair = pair1; // Assuming the components are copied; doesn't matter either way
            higherPair = pair2;
        }
        else
        {
            lowerPair = pair2;
            higherPair = pair1;
        }
        return higherPair.min - lowerPair.max;
    }

    // Compares the minima and maxima of two regions and returns whether they are overlapped or not
    // Returns true if overlapped; false otherwise
    bool AreRegionsBoundsOverlapped(Pair pair1, Pair pair2)
    {
        if(GetMinimumDistance(pair1, pair2) > OVERLAP_THRESHOLD) // Positive minimum distance: no overlap
            return false;
        else
            return true;
    }


    // Compares the minima and maxima of two regions and returns whether one region is nested by the other or not
    // Returns true is either region is nested in the other, and false otherwise
    bool ArePolymersNested(Pair pair1, Pair pair2)
    {
        Pair lowerPair, higherPair; // Pairs assigned based on the relative locations of the two regions

        // Assigning location:
        if(pair1.min < pair2.min)
        {
            lowerPair = pair1; // Assuming the components are copied; doesn't matter either way
            higherPair = pair2;
        }
        else
        {
            lowerPair = pair2;
            higherPair = pair1;
        }
        
        // Checking if the lower pair extends beyond the higher pair on teh other side as well:
        if(lowerPair.max > higherPair.max)
            return true;
        else
            return false;
    }
#endif