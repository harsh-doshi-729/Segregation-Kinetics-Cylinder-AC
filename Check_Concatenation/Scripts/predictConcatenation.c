// This script uses the maximal extension overlap method in order to predict whether the polymers and/or internal loops are concatenated or not

// TODO: Implement a method in region.h to get regions after reading a cross-links database.
// TODO: Use the above method to check pairwise concatenation of all the internal loops of an architecture.

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>

// importing paths file:
#include "../../Global_Scripts/System_File_Paths/system_file_paths.h"
#include GENERAL_FILE_METHODS
#include LAMMPS_POSITION_FILE
#include REGION
#include REGION_CONFIG

// Global variables:
int numberOfMonomers; // The number of monomers in a single polymer
int numberOfPolymers = 2; // The number of polymers in the system
int totalMonomers; // The total number of monomers in the system
char* architecture;

// Run information:
int runIndex = -1;

// Flag indicating whether internal concatenations should be checked for:
bool checkInternalConcatenations = false; // No by default

// File information:
char* directory; // The path to the directory where all the concatenation dump files are stored
simulation_data read_data; // The struct for the LAMMPS Position dump file

// region data:
region_data regionData; // The struct to store and read the information about regions in a polymer

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

// Accepts the arguments from the command line and sets them to global variables
void SetConstants(int argc, char** argv)
{
    int numberOfMandatoryArguments = 3;
    if(argc < numberOfMandatoryArguments + 1)
    {
        printf("Not enough arguments passed! Please pass the number of monomers (integer), the architecture(string), and the run index (integer) as arguments from the command line\n");
        printf("An optional flag (0 or 1) argument can be passed to indicate whether the internal concatenations should be checked too.\n");
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

        // Setting absolute path for directory:
        char* checkConcatenationPath;
        SetAbsolutePath(&checkConcatenationPath, CHECK_CONCATENATION);
        // Setting directory:
        SetDirectoryPath(&directory, checkConcatenationPath, numberOfMonomers);
        free(checkConcatenationPath);

        // Optional argument
        if(argc > numberOfMandatoryArguments + 1)
        {
            if(strcmp(argv[numberOfMandatoryArguments+1], "1") == 0)
                checkInternalConcatenations = true;
        }
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

// Sets the read file path to the dump file for the internal loops of a single polymer;
// The index of polymer starts from 1
void SetInternalLoopsReadFilePath(char** filePathPointer, int indexOfPolymer)
{
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    // Setting the file name according to the polymer index
    char* fileName;
    int bytes = asprintf(&fileName, "visual%i.dump", indexOfPolymer);
    if(bytes == -1)
    {
        printf("Memory could not be allocated for the file name! Terminating.\n");
        exit(1);
    }
    Append(filePathPointer, directory, folder, fileName);
    free(fileName);
    free(folder);
}

// Sets the values of the pair of min and max values for each region for each spatial dimension from the current timestep
void SetMinMax(Pair pairs[][3])
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

// Compares the minima and maxima of two regions and returns whether they are overlapped or not
// Returns true if overlapped; false otherwise
bool AreRegionsBoundsOverlapped(Pair pair1, Pair pair2)
{
    Pair lowerPair, higherPair; // Pairs assigned based on the relative locations of the two regions

    // Assigning location:
    if(pair1.min < pair2. min)
    {
        lowerPair = pair1; // Assuming the components are copied; doesn't matter either way
        higherPair = pair2;
    }
    else
    {
        lowerPair = pair2;
        higherPair = pair1;
    }

    double const OVERLAP_THRESHOLD = -1; // Threshold set to -sigma; there must be atleast a monomer distance between regions 
    // for them to be deemed unoverlapped
    // Comparing overlap parameter low.max-high.min with threshold:
    // overlap parameter = low.max - high.min; positive when there is overlap; otherwise negative
    if(lowerPair.max - higherPair.min > OVERLAP_THRESHOLD) // Positive overlap parameter, overlapped state
        return true;
    else
        return false;
}

// Checks whether the regions with the IDs are overlapped along all three dimensions; based on the pairs passed
// Returns true if they are, false if they are not overlapped along atleast one direction
// The regionIDs are indexed from 0
bool AreRegionsOverlapped(Pair pairs[][3], int regionID_1, int regionID_2)
{
    bool areOverlapped = true;
    for(int j = 0; j < 3; j++)
        areOverlapped = areOverlapped && AreRegionsBoundsOverlapped(pairs[regionID_1][j], pairs[regionID_2][j]);
    return areOverlapped; // true only if the regions are overlapped along all three directions
}

// Predicts the interpolymer concatenation by checking for overlap between the polymers at all time steps
bool PredictPolymerConcatenation(void)
{
    // Initializing the simulation_data struct:
    char* readFilePath;
    SetReadFilePath(&readFilePath);
    InitializeSimulationData(&read_data, totalMonomers, readFilePath);

    // Initializing region data:
    InitializeRegionManually(&regionData, numberOfPolymers, numberOfMonomers, true, false);
    // irrespective of value in regions_config file; so that polymers are separate regions

    // Reading the position data one timestep at a time:
    Pair pairs[regionData.totalNumberOfRegions][3];
    bool areConcatenated = true; // The flag that will be updated if the polymers stop being overlapped
    while(ReadPositions(&read_data))
    {
        SetMinMax(pairs);
        areConcatenated = AreRegionsOverlapped(pairs, 0, 1); // first and second regions: polymer regions
        if(!areConcatenated) // regions are no longer overlapped, hence unconcatenated
        {
            break;
        }
    }
    printf("%s Run %i: ", architecture, runIndex);
    if(areConcatenated)
        printf("Concatenated* Polymers!\n");
    else
        printf("Polymers Unconcatenated at iteration %i\n", read_data.timeStep);

    // Freeing up memory:
    FreeAllocatedMemory(&read_data);
    FreeRegionMemory(&regionData);
}

// Initializes the square matrix of flags to true; only the upper triangular part of the matrix is initialized
void InitializeFlagMatrix(size_t length, bool matrix[length][length])
{
    for(int i = 0; i < length; i++)
    {
        for(int j = i + 1; j < length; j++)
            matrix[i][j] = true;
    }
}

// Prints the upper diagonal elements of the square flag matrix; for debugging
void PrintUpperDiagonalMatrix(size_t length, bool matrix[length][length])
{
    for(int i = 0; i < length; i++)
    {
        for(int j = 0; j < length; j++)
        {
            if(j > i && i >= NON_LOOP_REGIONS)
                printf("%i ", matrix[i][j]);
            else
                printf("x ");
        }
        printf("\n");
    }
}

// Predicts the concatenation between different internal loops of a single polymer;
// The polymers are assumed to be indexed starting from 1
// TODO: Make it work for more complicated architectures with multiple disconnected loops
bool PredictInternalConcatenation(int indexOfPolymer)
{
    // Initializing the region data: assuming the regions are correctly set in the config
    InitializeRegion(&regionData, 1, numberOfMonomers); // looking at one polymer at a time

    // Initializing the simulation_data struct:
    char* readFilePath;
    SetInternalLoopsReadFilePath(&readFilePath, indexOfPolymer);
    // Initializing the struct as if all the monomers in a polymer are present
    InitializeSimulationData(&read_data, numberOfMonomers, readFilePath);
    // Suppressing warnings:
    read_data.suppressWarning = true;


    // The first few regions are assumed to be the polymer excluding internal loops which are deleted in these simulations
    // Setting up the 2D array of pairs for each pair of regions in three directions:
    int numberOfRegions = regionData.totalNumberOfRegions - NON_LOOP_REGIONS; // The actual number of regions excluding the non-looped regions
    int numberOfPairs = numberOfRegions * (numberOfRegions - 1) / 2; // N choose 2; for a particular direction
    Pair pairs[regionData.totalNumberOfRegions][3];

    // Checking for the pairwise concatenations between each pair of regions except the first one;
    bool areRegionsConcatenated[regionData.totalNumberOfRegions][regionData.totalNumberOfRegions]; // matrix of flags of concatenation between two regions
    InitializeFlagMatrix(regionData.totalNumberOfRegions, areRegionsConcatenated);
    bool areAllRegionsUnconcatenated; //  A flag to check if all regions are unoverlapped at a timestep
    while(ReadPositions(&read_data))
    {
        SetMinMax(pairs);

        areAllRegionsUnconcatenated = true; // resetting to default value
        // Iterating over all pairs of regions except for the non-looped ones:
        for(int i = NON_LOOP_REGIONS; i < regionData.totalNumberOfRegions; i++)
        {
            for(int j = i + 1; j < regionData.totalNumberOfRegions; j++)
            {
                if(areRegionsConcatenated[i][j]) // only checking the regions which were concatenated previously
                {
                    areRegionsConcatenated[i][j] = AreRegionsOverlapped(pairs, i, j);
                    if(areRegionsConcatenated[i][j]) // if the regions are overlapped
                        areAllRegionsUnconcatenated = false;
                }
            }
        }
        if(areAllRegionsUnconcatenated) // if none of the regions were overlapped
            break;
    }

    // Printing message about concatenations:
    printf("%s Run %i Polymer %i: ", architecture, runIndex, indexOfPolymer); // Printing prefix information for current case
    if(areAllRegionsUnconcatenated)
        printf("All internal regions unconcatenated!\n");
    else
        printf("At least one pair of internal regions was concatenated!\n");

    // Freeing up memory:
    FreeAllocatedMemory(&read_data);
    FreeRegionMemory(&regionData);

}

int test_main(int argc, char** argv)
{
    // Testing pass by reference behaviour of a struct
    Pair pairs[2][3];
    double array[][3] = {{1, 3, 2}, {7, 4, 5}, {3, 5, 3}, {9, 0, 7}, {2, 3, 1}};
    // reg1 : (1, 3), (3, 5), (1, 3)
    // reg2 : (7, 9), (0, 4), (5, 7)
    int size = 5;
    double** lookUpTable = malloc(size * sizeof(double*));
    for(int l = 0; l < size; l++)
        lookUpTable[l] = array[l];
    size_t length = 5;
    printf("Before Setting:\n");
    char directions[] = {'x', 'y', 'z'};
    for(int n = 0; n < 2; n++)
    {
        for(int j = 0; j < 3; j++)
        {
            printf("Region %i %c: \n", n + 1, directions[j]);
            PrintPair(pairs[n][j]);
        }
    }

    SetMinMax(pairs);

    printf("After Setting:\n");
    for(int n = 0; n < 2; n++)
    {
        for(int j = 0; j < 3; j++)
        {
            printf("Region %i %c: \n", n + 1, directions[j]);
            PrintPair(pairs[n][j]);
        }
    }
    free(lookUpTable);
    return 0;
}

int main(int argc, char** argv)
{
    SetConstants(argc, argv);
    // printf("Inter Polymer Concatenations:\n");
    PredictPolymerConcatenation();
    if(checkInternalConcatenations)
    {
        // printf("\nInternal Concatenations:\n");
        for(int n = 1; n <= numberOfPolymers; n++)
            PredictInternalConcatenation(n);
    }
    free(directory);
}