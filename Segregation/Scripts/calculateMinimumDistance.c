// This script predicts whether two polymers mixed in an infinite cylinder are segregated or not
// If they achieve segregation, the time of segregation is reported.
// It uses a very similar algorithm as used to predict concatenations between two polymers.

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>

// importing paths file:
#include "../../../../Scripts/System_File_Paths/system_file_paths.h"
#include GENERAL_FILE_METHODS
#include LAMMPS_POSITION_FILE
#include REGION
#include REGION_CONFIG
#include EXTENT

// Global variables:
int numberOfMonomers; // The number of monomers in a single polymer
int numberOfPolymers = 2; // The number of polymers in the system
int totalMonomers; // The total number of monomers in the system
char* architecture;
int runIndex;

// File information:
char* directory; // The path to the directory where all the concatenation dump files are stored
simulation_data read_data; // The struct for the LAMMPS Position dump file

// region data:
region_data regionData; // The struct to store and read the information about regions in a polymer

// Flag to indicate whether to use regions defined in the region_config.h file
bool useRegions = false; // Default value
int regionIDs[2] = {1, 2}; // The IDs (indexed from 1) for the regions between which the minimum distance is to be calculated

// Accepts the arguments from the command line and sets them to global variables
void SetConstants(int argc, char** argv)
{
    int numberOfMandatoryArguments = 3;
    if(argc < numberOfMandatoryArguments + 1)
    {
        printf("Not enough arguments passed! Please pass the number of monomers (integer), the architecture(string), and the run index (integer) as arguments from the command line\n");
        printf("For example: ./a.out 200 Arc1_10 1\n");
        printf("An optional boolean argument can be passed after the above arguments to indicate whether the regions in regions_config.h should be used\n");
        printf("If the optional argument passed is true/1, two region IDs should also be passed after the flag argument.\n");
        exit(1);
    }
    else
    {
        numberOfMonomers = atoi(argv[1]);
        if(numberOfMonomers == 0)
        {
            printf("The argument passed (%s) as number of monomers could not be converted to a valid number! Terminating.\n", argv[1]);
            exit(1);
        }
        totalMonomers = numberOfMonomers * numberOfPolymers;

        architecture = argv[2];

        runIndex = atoi(argv[3]);
        if(runIndex == 0)
        {
            printf("The argument passed (%s) as the run index could not be converted to a valid number! Terminating.\n", argv[3]);
            exit(1);
        }

        // Setting directory:
        SetDirectoryPath(&directory, NEW_SEGREGATION, numberOfMonomers);

        // Getting the optional argument if passed:
        if(argc > numberOfMandatoryArguments + 1)
        {
            char* regionFlag = argv[numberOfMandatoryArguments+1];
            if(strcasecmp(regionFlag, "true") == 0 || strcmp(regionFlag, "1") == 0)
                useRegions = true;
            else if(strcasecmp(regionFlag, "false") == 0 || strcmp(regionFlag, "0") == 0)
                useRegions = false;
            else
            {
                printf("The optional argument passed (%s) to indicate whether to use regions is not recognised! Please pass either true/1 or false/0\n", regionFlag);
                exit(1);
            }
            // Reading region IDs:
            if(useRegions)
            {
                if (argc < numberOfMandatoryArguments + 3)
                {
                    printf("ERROR: The region IDs were not passed after the optional flag! Please pass two integers for the region IDs between which the minimum distance can be calculated.\n");
                    exit(1);
                }
                else
                {
                    regionIDs[0] = atoi(argv[numberOfMandatoryArguments+2]);
                    regionIDs[1] = atoi(argv[numberOfMandatoryArguments+3]);
                    if(regionIDs[0] == 0 || regionIDs[1] == 0)
                    {
                        printf("ERROR: At least one of the passed arguments for region IDs ('%s', '%s') is invalid. Please pass positive integers for the region IDs.\n", argv[numberOfMandatoryArguments+2], argv[numberOfMandatoryArguments+3]);
                        exit(1);
                    }
                }
            }
        }
    }
}

// Verifies the validity of the direction and region indices
void VerifyIndices(int directionIndex, int regionID1, int regionID2)
{
    if(directionIndex < 0 || directionIndex > 2) // Direction index must be 0 (x), 1 (y), or 2 (z)
    {
        printf("ERROR: AreRegionsOverlapped(): The directionIndex %i is out of bounds.\n", directionIndex);
        exit(1);
    }
    if(regionID1 < 0 || regionID2 < 0 || regionID1 > regionData.totalNumberOfRegions-1 || regionID2 > regionData.totalNumberOfRegions-1)
    {
        printf("ERROR: AreRegionsOverlapped(): At least one of the regionsIDs (%i, %i) is out of bounds. Number of regions: %i.\n", regionID1, regionID2, regionData.totalNumberOfRegions);
        exit(1);
    }
}

// Sets the read file path of the dump file
void SetReadFilePath(char** filePathPointer)
{
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    Append(filePathPointer, directory, folder, "visual.dump");
    free(folder);
}

// Sets the file path to write the minimum distance trajectory
void SetMinDistanceFilePath(char** filePathPointer)
{
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    char* fileName = "min_distance.dat";
    if(useRegions)
    {
        int bytes = asprintf(&fileName, "min_distance_reg_%i_%i.dat", regionIDs[0], regionIDs[1]);
        if(bytes == -1)
        {
            printf("ERROR: Memory could not be allocated for the minimum distance file name with region IDs.\n");
            exit(1);
        }
    }
    Append(filePathPointer, directory, folder, fileName);
    free(folder);
    if(useRegions)
        free(fileName);
}

// Sets the file path to write the induction times for all runs
void SetInductionTimeFilePath(char** filePathPointer)
{
    char* fileName = "/Analysis/inductionTimes.csv";
    Append(filePathPointer, directory, architecture, fileName);
}

// Predicts whether the two regions are overlapped along one direction direction
bool AreRegionsOverlapped(Pair regionPairs[][3] , int directionIndex, int regionID1, int regionID2)
{
    // verifying validity of arguments:
    VerifyIndices(directionIndex, regionID1, regionID2);

    // Pairs to compare:
    Pair pair1 = regionPairs[regionID1][directionIndex];
    Pair pair2 = regionPairs[regionID2][directionIndex];
    return AreRegionsBoundsOverlapped(pair1, pair2);
}

// Predicts whether the polymers segregate and returns the time of segregation.
// If the polymers do not segregate, returns -1
void PredictSegregation(void)
{
    // Initializing the simulation_data struct:
    char* readFilePath;
    SetReadFilePath(&readFilePath);
    InitializeSimulationData(&read_data, totalMonomers, readFilePath);

    // Initializing region data:
    if(!useRegions) // ignoring regions_config file; treating each polymer as a separate region
        InitializeRegionManually(&regionData, numberOfPolymers, numberOfMonomers, true, false);
    else
    {
        printf("Reading regions from regions_config file. Please ensure the regions are set as desired.\n");
        InitializeRegion(&regionData, numberOfPolymers, numberOfMonomers);
    }
    bool arePolymersOverlapped = true; // The polymers are expected to be overlapped initially
    // TODO: Compute whether polymers are nested, and calculate induction time

    // Reading the position data one timestep at a time:
    Pair pairs[regionData.totalNumberOfRegions][3];
    while(ReadPositions(&read_data))
    {
        SetMinMax(pairs, regionData, read_data);
        arePolymersOverlapped = AreRegionsOverlapped(pairs, 2, regionIDs[0] - 1, regionIDs[1] - 1); // first and second regions: polymer regions
        if(!arePolymersOverlapped) // regions are no longer overlapped, have segregated
        {
            break;
        }
    }
    printf("%s Run %i: ", architecture, runIndex);
    if(arePolymersOverlapped)
        printf("Polymers not segregated!\n");
    else
        printf("Polymers Segregated at iteration %i\n", read_data.timeStep);

    // Freeing up memory:
    FreeAllocatedMemory(&read_data);
    FreeRegionMemory(&regionData);
}

// Writes the minimum distance between two regions for the current timestep to the passed file
void WriteMinimumDistance(FILE* writeFilePointer, Pair pairs[][3], int directionIndex, int regionID1, int regionID2)
{
    // verifying validity of arguments:
    VerifyIndices(directionIndex, regionID1, regionID2);
    // Pairs to compare:
    Pair pair1 = pairs[regionID1][directionIndex];
    Pair pair2 = pairs[regionID2][directionIndex];

    double minDistance = GetMinimumDistance(pair1, pair2);
    bool areRegionsOverlapped = AreRegionsBoundsOverlapped(pair1, pair2);

    // Writing to file:
    fprintf(writeFilePointer, "%i, %lf, %i\n", read_data.timeStep, minDistance, areRegionsOverlapped);
}

// Reads the positions and writes the minimum distance trajectory along for each timestep
// Also calculates and writes the induction time of the two polymers
void IterateOverAllTimesteps(void)
{
    // Initializing the simulation_data struct:
    char* readFilePath;
    SetReadFilePath(&readFilePath);
    InitializeSimulationData(&read_data, totalMonomers, readFilePath);

    // Initializing region data:
    if(!useRegions) // ignoring regions_config file; treating each polymer as a separate region
        InitializeRegionManually(&regionData, numberOfPolymers, numberOfMonomers, true, false);
    else
    {
        printf("Reading regions from regions_config file. Please ensure the regions are set as desired.\n");
        InitializeRegion(&regionData, numberOfPolymers, numberOfMonomers);
    }

    // Setting up file to write minimum distance:
    char* writeFilePath;
    SetMinDistanceFilePath(&writeFilePath);
    FILE* writeFilePointer = fopen(writeFilePath, "w");

    bool arePolymersOverlapped = true; // The polymers are expected to be overlapped initially

    // Defining variables for calculating induction time:
    bool arePolymersNested = false; // Assuming the polymers are not nested; default
    bool isNestedStateInitialized = false; // A flag to indicate whether teh arePolymersNested flag has been initialized based on the first snapshot of the run
    int inductionTime = -1; // the time it takes for the polymers to go from a nested state to un-nested state
    // Assumption: The polymers can only be nested if they start out nested. 
    // It is entropically unfavourable for the polymers to become nested from an un-nested state
    // However, teh latter can occur as fluctuations; which we disregard

    // Writing file header:
    fprintf(writeFilePointer, "Timestep, minimum distance, are overlapped\n");

    // Reading the position data one timestep at a time:
    Pair pairs[regionData.totalNumberOfRegions][3];
    while(ReadPositions(&read_data))
    {
        // Setting the mix-max pair values for current timestep:
        SetMinMax(pairs, regionData, read_data);
        // Initializing nested state if not already done
        if(!isNestedStateInitialized)
        {
            arePolymersNested = ArePolymersNested(pairs[regionIDs[0]-1][2], pairs[regionIDs[1]-1][2]); // Checking nestedness between polymers along z axis
            isNestedStateInitialized = true;
        }
        else if(arePolymersNested) // If the polymers were nested; checking if they became un-nested this timestep
        {
            arePolymersNested = ArePolymersNested(pairs[regionIDs[0]-1][2], pairs[regionIDs[1]-1][2]);
            if(!arePolymersNested) // The polymers stopped being nested
                inductionTime = read_data.timeStep;
        }
        WriteMinimumDistance(writeFilePointer, pairs, 2, regionIDs[0] - 1, regionIDs[1] - 1); // Writing minimum distance along z direction between the two regions (polymer regions)
    }

    // Printing induction time information:
    if(!useRegions) // Only printing induction time if polymer regions are used
    {
        if(inductionTime == -1) // Induction period absent
            printf("The induction period was absent for %s Run %i. The polymers were not nested.\n", architecture, runIndex);
        else
        {
            printf("Induction time for %s Run %i: %i\n", architecture, runIndex, inductionTime);
            // Writing to file:
            char* inductionFilePath;
            SetInductionTimeFilePath(&inductionFilePath);
            FILE* inductionFilePointer = fopen(inductionFilePath, "a"); // Appending to the file
            fprintf(inductionFilePointer, "%i, %i\n", runIndex, inductionTime);
            fclose(inductionFilePointer);
            free(inductionFilePath);
        }
    }

    // Closing write file:
    fclose(writeFilePointer);
    // Freeing up memory:
    free(writeFilePath);
    FreeAllocatedMemory(&read_data);
    FreeRegionMemory(&regionData);
}

// Predict Segregation for all runs

int main(int argc, char** argv)
{
    SetConstants(argc, argv);
    // PredictSegregation();
    IterateOverAllTimesteps();
}