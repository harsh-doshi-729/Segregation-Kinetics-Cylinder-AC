// This header file contains variables and functions to help divde aa polymer into different regions.
// These regions are usually different loops in the same polymer, but can be set to anything else

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifndef REGIONS_H // include guard
    #define REGIONS_H
    // Importing the region bounds from a config file:
    #include "../Config_Files/regions_config.h"
    #include "../System_File_Paths/system_file_paths.h"
    #include MODULO

    // Defining a struct for storing the upper and lower bounds (monomer indices) of a region
    struct region_bounds {
        int upper;
        int lower;
    };

    typedef struct region_bounds region_bounds;

    // Defining a struct to store all the region related information:
    struct region_data {
        int numberOfPolymers; // The number of polymers in the system
        int numberOfMonomers; // The number of monomers in a single polymer
        int numberOfRegions; // The number of regions in a single polymer
        int totalNumberOfRegions; // The total number of regions over all polymers in the system

        region_bounds* regions; // An array of regions bounds demarcating the regions for the first polymer

        bool doNotUseRegions; // A flag to indicate whether to divide polymers into subregions. 
        // If true, the polymers are divided into subregions according to the regions array. If false, each polymer is treated as a separate region

        bool divideOnlyFirstPolymer; // If true, only the first polymer is divided into subregions according to the regions array 
        // while the other polymers are treated as individual regions without subdivision. If false, all polymers are subdivided into regions
        // This is useful when the first polymer needs to be analysed more closely, but not the others. For example: in the Loop Cutting problem
    };

    typedef struct region_data region_data;

    // Functions to set up the struct:

    // Allocates memory for the region_bounds array
    void AllocateRegionBounds(region_data* data, size_t arrayLength)
    {
        data->regions = (region_bounds*) malloc(arrayLength * sizeof(region_bounds));
        if (data->regions == NULL)
        {
            printf("Memory could not be allocated for the region bounds! Terminating.\n");
            exit(1);
        }
    }

    // Sets the default regions as each polymer. Used when doNotUseRegions flag is true
    void SetDefaultRegions(region_data* data)
    {
        AllocateRegionBounds(data, data->totalNumberOfRegions);
        for(int i = 0; i < data->totalNumberOfRegions; i++)
        {
            data->regions[i].lower = i * data->numberOfMonomers + 1; // The first monomer in the polymer
            data->regions[i].upper = (i+1) * data->numberOfMonomers; // the last monomer in the polymer
        }
    }

    // TODO: Add a way to set region bounds based on cross links

    // Allocates memory and initializes the array of region bounds by reading the regions from the config file. Returns the number of regions read.
    int AllocateAndSetRegionBounds(region_data* data)
    {
        int regionArray[REGIONS_ARRAY_LENGTH][2] = REGION_BOUNDS_ARRAY; // Defining the array from the values listed in the config file
        AllocateRegionBounds(data, REGIONS_ARRAY_LENGTH);
        for(int i = 0; i < REGIONS_ARRAY_LENGTH; i++)
        {
            data->regions[i].lower = regionArray[i][0];
            data->regions[i].upper = regionArray[i][1];
        }
        return REGIONS_ARRAY_LENGTH;
    }

    // Sets the values of region bounds and other variables after the necessary variables have been set
    void SetRegion(region_data* data)
    {
        // Setting the values of the other member variables:
        if(data->doNotUseRegions)
        {
            data->numberOfRegions = 1; // 1 region per polymer
            data->totalNumberOfRegions = data->numberOfPolymers;
            SetDefaultRegions(data);
        }
        else
        {
            data->numberOfRegions = AllocateAndSetRegionBounds(data); // The number of regions in the first polymer
            data->totalNumberOfRegions = data->numberOfRegions; // Initializing total regions with the value of the first polymer; updated later
            if(data->divideOnlyFirstPolymer)
            {
                data->totalNumberOfRegions += data->numberOfPolymers - 1; // Adding one region for each of the remaining polymers
            }
            else
            {
                data->totalNumberOfRegions += (data->numberOfPolymers - 1) * data->numberOfRegions; // Adding the same number of regions for each of the remaining polymers
            }
        }
    }

    // Initializes the region_data struct with the values passed and the regions read from a config file
    void InitializeRegion(region_data* data, int numberOfPolymers, int numberOfMonomers)
    {   
        data->numberOfPolymers = numberOfPolymers;
        data->numberOfMonomers = numberOfMonomers;
        data->doNotUseRegions = DO_NOT_USE_REGIONS;
        data->divideOnlyFirstPolymer = DIVIDE_ONLY_FIRST_POLYMER; // taking the flag values from the config file

        SetRegion(data);
    }

    // Initializes the region_data struct manually from the passed values (not from the config file)
    void InitializeRegionManually(region_data* data, int numberOfPolymers, int numberOfMonomers, bool doNotUseRegions, bool divideOnlyFirstPolymer)
    {
        data->numberOfPolymers = numberOfPolymers;
        data->numberOfMonomers = numberOfMonomers;
        data->doNotUseRegions = doNotUseRegions;
        data->divideOnlyFirstPolymer = divideOnlyFirstPolymer;

        SetRegion(data);
    }

    // Frees up the memory allocated for the region struct:
    void FreeRegionMemory(region_data* data)
    {
        free(data->regions);
    }

    // Utility Functions:

    // Returns the region ID according to the region definitions by the regions array; indexing starts from 0
    // The monomer index is expected to start from 1
    int GetRegionID(region_data data, int monomerIndex) // region_data argument is read-only
    {
        int polymerIndex = (monomerIndex - 1) / data.numberOfMonomers; // the polymer index with base 0
        int reducedMonomerID = (monomerIndex - 1) % data.numberOfMonomers + 1; // The index of the monomer in the first polymer corresponding to the passed monomer

        if(data.doNotUseRegions)
            return polymerIndex; // returning the polymer index as the region in case the read regions do not need to be used

        int regionID = -1; // default
        for(int i = 0; i < data.numberOfRegions; i++)
        {
            if (data.regions[i].lower <= data.regions[i].upper) // normal case
            {
                // the monomer lies directly between the lower and upper bounds
                if(reducedMonomerID >= data.regions[i].lower && reducedMonomerID <= data.regions[i].upper)
                {
                    regionID = i;
                    break;
                }
            }
            else // abnormal case, the discontinuity in monomer index occurs in this region
            {
                // the monomer lies between 1 and upper or between lower and numberofMonomers
                if(reducedMonomerID >= data.regions[i].lower || reducedMonomerID <= data.regions[i].upper)
                {
                    regionID = i;
                    break;
                }
            }
        }

        // Checking if the region was found:
        if(regionID == -1)
        {
            printf("The region ID for monomer %i could not be found! Please check the region bounds in the config file. Terminating.\n", monomerIndex);
            exit(1);
        }

        if(data.divideOnlyFirstPolymer && polymerIndex > 0)
            return data.numberOfRegions + polymerIndex - 1; // -1 since the regionID is 0 based
        else // if all polymers are divided or the monomer belongs to the first polymer
            return data.numberOfRegions * polymerIndex + regionID; // The region in that particular polymer
    }

    // Returns the number of monomers in one of the regions demarcated by region_bounds
    // The regionIndex is expected to be indexed starting from 0
    int GetRegionCount(region_data data, int regionIndex) // region_data is read-only (copy only)
    {
        int polymerIndex = regionIndex / data.numberOfRegions;
        regionIndex = modulo(regionIndex, data.numberOfRegions); // converting the region index to the corresponding region index in the first polymer
        // Checking if the region is the entire polymer:
        if(data.doNotUseRegions || data.divideOnlyFirstPolymer && polymerIndex > 0)
            return data.numberOfMonomers; // All monomers in the polymer belong to one region
        else
            return modulo(data.regions[regionIndex].upper - data.regions[regionIndex].lower, data.numberOfMonomers) + 1;
            // The modulo makes it so that this works even if the upper and lower indices are separated by the ori discontinuity i.e the 200-1 bond.
    }

    // Returns the number of monomers present in a particular region in one of the polymers
    // The regionID is expected to be indexed starting from 0
    int GetNumberOfMonomersInRegion(region_data data, int regionID)
    {
        // Assuming the regionID passed is based 0
        int polymerIndex = regionID / data.numberOfRegions;
        int reducedRegionID = regionID % data.numberOfRegions;

        if(data.doNotUseRegions || data.divideOnlyFirstPolymer && polymerIndex > 0)
            return data.numberOfMonomers; // The polymer is an entire region
        else
        {
            return GetRegionCount(data, reducedRegionID);
        }
    }

    // Accepts a simulation_data struct to calculate the CoM of the monomers belonging to one region specified by regionID (indexed from 0)
    // The three CoM coordinates are calculated and stored in the passed array
    void CalculateRegionCoM(double* CoM, simulation_data read_data, region_data regData, int regionID)
    {
        int counter = 0;
        // Initializing:
        for(int j = 0; j < 3; j++)
        {
            CoM[j] = 0;
        }
        // Calculating CoMs
        for(int i = 0; i < read_data.numberOfMonomers; i++)
        {
            int queryRegion = GetRegionID(regData, i+1);
            if(queryRegion == regionID)
            {
                for(int j = 0; j < 3; j++)
                    CoM[j] += read_data.bead_positions[i][j];
                counter++;
            }
        }
        // Normalizing:
        if(counter != GetRegionCount(regData, regionID))
        {
            printf("ERROR: The region count does not match the number of monomers counted while calculating region CoM fo region %i.\n", regionID + 1);
            exit(1);
        }
        // else:
        for(int j = 0 ; j < 3; j++)
        {
            CoM[j] /= counter;
        }
    }
#endif