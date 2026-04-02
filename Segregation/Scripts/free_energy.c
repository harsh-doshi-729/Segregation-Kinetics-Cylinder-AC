// This script calculates the free energy of different configurations of segregated polymers with respect to a reaction coordinate.
// The free energy is back-calculated from the probability of the system to exist with a value of the reaction coordinate.


#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

// importing paths file:
#include "../../../../Scripts/System_File_Paths/system_file_paths.h"
#include GENERAL_FILE_METHODS
#include LAMMPS_POSITION_FILE
#include MONOMER_DISTRIBUTION // To calculate the probability distribution with respect to the reaction coordinate
#include REGION
#include SEGREGATION_TIMES // To read the segregation times file
#include CROSS_LINKS_DATABASE // To read the cylinder diameter 
#include EXTENT // To calculate overlap distance

// Steps:
// 1. Read the data: The position files of the segregation runs must be read such that only data after segregation is used. Count number of steps
// 2. Calculate probability distribution with respect to a reaction coordinate
// 3. Invert to calculate free energy.
// 4. Set appropriate offset/constant for the free energy

// Global Variables:
int numberOfPolymers = 2; // The total number of polymers in the system
int numberOfMonomers; // The number of monomers in a single polymer of the system
char* architecture; // The name of the architecture being used
int numberOfPairs; // The number of pairs of regionIDs; each pair defines a reaction coordinate: distance between CoMs of the regions
int** regionPairs; // 2D Array to store the pairs of regions between whom CoM distance is to be calculated as the reaction coordinate
const int direction = 2; // z direction

// Directory information:
char* directory;
int runIndex = -1;

// Segregation Times file information:
int numberOfRuns = 50;
segregation_times segData;

// Position file struct:
simulation_data read_data;

// Struct to calculate the probability distribution with respect to the reaction coordinate
monomer_distribution_data dist_data;
int numberOfBins; // The number of bins for the distribution function; accepted from the user
double lowerBoundary; // The lower limit for the reaction coordinate;
double upperBoundary; // The upper limit for the reaction coordinate; set in accordance with the reaction coordinate chosen

// Flag to indicate whether the unsegregated data should be skipped while calculating the probability distribution
bool skipUnsegregatedData = false;
// Flag to indicate whether the free energy should be calculated and printed:
bool calculateFreeEnergy = false; // should only be calculated when the simulation is an equilibrium simulation

// Region struct to distnguish monomers into different regions:
region_data regData;

// Struct to store the cylinder dimensions for this architecture:
Architecture archDimensions;

// Mode of operation: the reaction coordinate with respect to which the distribution and/or free energy should be calculated
char* acceptedReactionCoords[] = {"region_CoM", "region_overlap"};
int numberOfAcceptedReactCoords = 2;
int chosenReactionCoord;

// Allocating memory for the region pairs:
void AllocateRegionPairs(int numberOfPairs)
{
    regionPairs = (int**) malloc(numberOfPairs * sizeof(int[2]));
    for(int i = 0; i < numberOfPairs; i++)
    {
        regionPairs[i] = (int*) malloc(2 * sizeof(int));
    }
}

// Frees the memory allocated to regionPairs
void FreeRegionPairs(void)
{
    for(int i = 0 ; i < numberOfPairs; i++)
        free(regionPairs[i]);
    free(regionPairs);
}

void SetConstants(int argc, char** argv)
{
    int numberOfMandatoryArguments = 6;
    if(argc < numberOfMandatoryArguments + 1)
    {
        printf("Not enough arguments passed! Please pass the number of monomers (integer), the architecture name (string), the number of bins for the free energy or distribution (integer), a flag to indicate whether the free energy should be calculated, the reaction coordinate identifier, the number of regionID pairs (integer) followed by the corresponding pairs of regionIDs (space-separated integers) as arguments to the command line.\n");
        printf("Accepted reaction coordinate identifiers: ");
        PrintArray(numberOfAcceptedReactCoords, acceptedReactionCoords);
        printf("For example: ./a.out 200 Arc2 100 true region_CoM 1 1 2\n");
        printf("An optional argument for the run index can be passed after these arguments to perform the calculation for only that run.\n");
        exit(1);
    }
    // else:
    // Number of monomers:
    numberOfMonomers = atoi(argv[1]);
    if(numberOfMonomers == 0)
    {
        printf("The argument %s could not be converted to a valid integer! Terminating.\n", argv[1]);
        exit(1);
    }


    // Architecture
    architecture = argv[2];

    // Number of bins:
    numberOfBins = atoi(argv[3]);
    if(numberOfBins == 0)
    {
        printf("The argument %s could not be converted to a valid integer! Terminating.\n", argv[3]);
        exit(1);
    }
    
    // Setting flag:
    if(argv[4] == "true" || atoi(argv[4]) == 1)
    {
        skipUnsegregatedData = true;
        calculateFreeEnergy = true;
    }
    else if(argv[4] == "false" || atoi(argv[4]) == 0)
    {
        skipUnsegregatedData = false;
        calculateFreeEnergy = false;
    }

    char* reactionCoordArg = argv[5];
    chosenReactionCoord = SearchStringArray(reactionCoordArg, numberOfAcceptedReactCoords, acceptedReactionCoords);
    if(chosenReactionCoord == -1)
    {
        printf("ERROR: The reaction coordinate idnetifier '%s' was not recognized. Please pass one of the following accpeted indetifiers: ", reactionCoordArg);
        PrintArray(numberOfAcceptedReactCoords, acceptedReactionCoords);
        exit(1);
    }

    // Reading regionIDs:
    numberOfPairs = atoi(argv[6]);
    if(numberOfPairs <= 0)
    {
        printf("The argument %s for number of pairs could not be converted to a valid integer! Terminating.\n", argv[6]);
        exit(1);
    }
    int argumentsAccepted = numberOfMandatoryArguments; // The number of arguments already read and stored
    int numberOfIDs = 2 * numberOfPairs; // updating to include regionID pairs
    numberOfMandatoryArguments = numberOfMandatoryArguments + numberOfIDs; // updating number of mandatory arguments
    AllocateRegionPairs(numberOfPairs);
    for(int i = 0; i < numberOfIDs; i++)
    {
        int regionID = atoi(argv[argumentsAccepted + 1 + i]);
        if(regionID == 0)
        {
            printf("The argument %s for a regionID could not be converted to a valid integer! Terminating.\n", argv[argumentsAccepted + 1 + i]);
            exit(1);
        }
        int pairIndex = i / 2; // The index of the current pair
        int ID_index = i % 2; // The index of the regionID within a pair
        regionPairs[pairIndex][ID_index] = regionID;
    }

    // Upper Boundary of distributions: axisLength of the cylinder
    // Reading diameter and setting axis length
    InitializeArcConfinementDimensions(&archDimensions, numberOfMonomers, architecture);
    if(IsCylinderInfinite())
        upperBoundary = GetInfiniteCylinderAxisLength(numberOfMonomers);
    else
        upperBoundary = 5 * archDimensions.confinementDiameter;
    FreeArchitecture(&archDimensions);
    // Lower Boundary: 0 for COM distance
    lowerBoundary = 0; // Default
    // For Overlap distance, majority of the values will be negative; minimum value is negative axis length
    if(strcmp(acceptedReactionCoords[chosenReactionCoord], "region_overlap") == 0)
        lowerBoundary = -upperBoundary; // negative axis length
        upperBoundary = 2 * archDimensions.axisLength; // Assuming at full overlap, the polymers are stretched a maximum of twice of a single polymer

    // Optional argument:
    if(argc > numberOfMandatoryArguments + 1) // Atleast one optional argument passed
    {
        runIndex = atoi(argv[numberOfMandatoryArguments + 1]);
        if(runIndex == 0)
        {
            printf("The argument %s could not be converted to a valid integer! Terminating.\n", argv[numberOfMandatoryArguments + 1]);
            exit(1);
        }
    }
    // Setting the directory:
    SetDirectoryPath(&directory, NEW_SEGREGATION, numberOfMonomers);
}

// Sets the path to the position dump file
void SetReadFilePath(char** filePathPointer, int runIndex)
{
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    char* fileName = "visual.dump";
    Append(filePathPointer, directory, folder, fileName);
    free(folder);
}

// Sets the path to the file where the probability distribution of the reaction coordinate will be written
void SetDistributionFilePath(char** filePathPointer, int runIndex, char* label)
{
    char* fileNameSuffix = "_distribution.csv";
    // Handling the case of no run index passed:
    if(runIndex == -1)
    {
        runIndex = 1; // storing the file in the run1 folder
        fileNameSuffix = "_distribution_all.csv";
    }
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    char* fileName;
    Append(&fileName, label, fileNameSuffix, "");
    Append(filePathPointer, directory, folder, fileName);
    free(fileName);
    free(folder);
}

// Sets the free energy label (file name) for the free energy file incorporating the regions between the CoM distance is calculated
// The two region IDs accepted as arguments are expected to be indexed from 1
void SetFreeEnergyLabel(char** labelPointer, char* reactionCoordIdentifier, int region1, int region2)
{
    char* regionLabel; // A label introduced to avoid ambiguitiy between polymeric regions and sub-polymer regions
    if(DO_NOT_USE_REGIONS)
        regionLabel = "pol"; // The regions are polymers, or 'pol' for short
    else
        regionLabel = "reg"; // The usual regions that are a sub-contour of each polymer
    int bytes = asprintf(labelPointer, "%s_%s_%i_%i", reactionCoordIdentifier, regionLabel, region1, region2);
    if(bytes == -1)
    {
        printf("ERROR: Memory could not be allocated for the free energy label!\n");
        exit(1);
    }
}

// Sets the path to the file where the free energy as a function of the reaction coordinate will be written
// The label argument is to uniquely characterise the reaction coordinate with which the free energy is calculated
void SetFreeEnergyFilePath(char** filePathPointer, int runIndex, char* label)
{
    char* folder;
    if(runIndex != -1) // runIndex set
    {
        SetFolderName(&folder, architecture, runIndex);
    }
    else
    {
        Append(&folder, architecture, "/Analysis/", "FreeEnergy/"); // Saving it in the analysis folder
    }
    char* fileNamePrefix = "free_energy_";
    char* fileName;
    Append(&fileName, fileNamePrefix, label, ".csv");
    Append(filePathPointer, directory, folder, fileName);
    free(folder);
    free(fileName);
}

// Returns the distance between the CoMs of two regions in the system along the z axis
// Takes two integers as arguments for the two region IDs; indexed from 1
double GetRegionsZCoMDistance(simulation_data read_data, region_data regData, int region1, int region2) // Both arguments are read_only (only copies)
{
    int numberOfRegions = 2;
    int regionIDs[] = {region1, region2};
    double z_CoMs[numberOfRegions]; // An array to store the z coordinate of the CoMs of the two regions

    // Calculating CoM:
    for(int j = 0; j < numberOfRegions; j++)
    {
        double CoM[3];
        CalculateRegionCoM(CoM, read_data, regData, regionIDs[j] - 1);
        z_CoMs[j] = CoM[2]; // only storing the z coordinate
    }

    // returning z distance:
    return fabs(z_CoMs[0] - z_CoMs[1]);
}

// Returns the 'Minimum Distance' or overlap distance between two regions
// Args:
// - read_data: The simulation_data struct that contains the read monomer positions for the current timestep
// - regData: The region_data struct that stores information about the regions the system is divided in
// - region1: The index of the first region; indexed from 1
// - region2: The index of the second region; indexed from 1
double GetRegionsOverlapDistance(simulation_data read_data, region_data regData, int region1, int region2)
{
    Pair regionExtentPairs[regData.totalNumberOfRegions][3];
    SetMinMax(regionExtentPairs, regData, read_data);
    return -GetMinimumDistance(regionExtentPairs[region1 - 1][direction], regionExtentPairs[region2 - 1][direction]);
    // Overlap distance is the negative of the minimum distance
}

// Returns the value of the reaction coordinate for the timestep read in the passed simulation_data
// The reaction coordinate is chosen by the reaction coordinate identifier passed as an argument to this script
// Args:
// - read_data: The simulation_data struct that contains the read monomer positions for the current timestep
// - regData: The region_data struct that stores information about the regions the system is divided in
// - region1: The index of the first region; indexed from 1
// - region2: The index of the second region; indexed from 1
double GetReactionCoordinate(simulation_data read_data, region_data regData, int region1, int region2) // The read_data and regData are read-only
{
    if(strcmp(acceptedReactionCoords[chosenReactionCoord], "region_CoM") == 0)
        return GetRegionsZCoMDistance(read_data, regData, region1, region2);
    else if(strcmp(acceptedReactionCoords[chosenReactionCoord], "region_overlap") == 0)
        return GetRegionsOverlapDistance(read_data, regData, region1, region2);
}

// Prints the calculated probability distribution of the reaction coordinate to a file
// Args:
// - dist_data: A pointer to the monomer_distribution struct that contains the information of the probability distribution to be printed
// - distribution: A double array that contains the values of the calculated probability distribution
// - writeFilePrefix: The label to be appended to the write file path indicating the reaction coordinate for the distribution
void PrintProbabilityDistribution(monomer_distribution_data* dist_data, double* distribution, char* writeFilePrefix)
{
    // Printing probability distribution to file: debugging
    char* writeFilePath;
    SetDistributionFilePath(&writeFilePath, runIndex, writeFilePrefix);
    FILE* writeFile = fopen(writeFilePath, "w");
    if(writeFile == NULL)
    {
        printf("ERROR: The file %s could not be opened!\n", writeFilePath);
        exit(1);
    }
    free(writeFilePath);
    PrintDistribution(dist_data, distribution, true, writeFile);
    fclose(writeFile);
}

// Reads the position dump file, calculates and prints a histogram (probability distribution) of the reaction coordinate for the system
// The reaction coordinate is taken as the CoM difference between the two regions specified by the pair of region indices passed (indexed from 1)
// The distribution array passed must have been allocated memory
void ComputeProbabilityDistribution(size_t distributionLength, double distribution[distributionLength], int region1, int region2)
{
    // Initializing the monomer distribution data struct:
    InitializeMonomerDistributionDataWithNumberOfBins(&dist_data, numberOfMonomers, lowerBoundary, upperBoundary, numberOfBins);

    int runIterator; // A counter to iterate ove all runs from which the distribution is to be calculated
    int numberOfCalculationRuns; // The number of runs over which the distribution should be calculated
    if(runIndex == -1) // A specific runIndex has not been passed; use all runs
    {
        runIterator = 1; // The calculation will start with run number 1
        numberOfCalculationRuns = numberOfRuns; // Calculation will be performed for all runs
    }
    else
    {
        runIterator = runIndex;
        numberOfCalculationRuns = 1;
    }

    // Debugging: counting the number of snapshots that contribute to the distribution
    int counter = 0;
    do // iterating over atleast 1 run:
    {
        // Initializing the position dump file:
        char* readFilePath;
        SetReadFilePath(&readFilePath, runIterator);
        InitializeSimulationData(&read_data, numberOfPolymers * numberOfMonomers, readFilePath);
        if(!(skipUnsegregatedData && segData.segregationTimes[runIterator-1] == -1)) // segregation occured
        {
            // Reading the dump file and adding to the distribution:
            printf("Reading positions for Run %i.\n", runIterator);
            while(ReadPositions(&read_data))
            {
                // Skipping the timesteps before segregation:
                if((skipUnsegregatedData && read_data.timeStep <= segData.segregationTimes[runIterator-1]))
                    continue;
                
                IncrementMonomerCounts(&dist_data, GetReactionCoordinate(read_data, regData, region1, region2));
                counter++;
            }
        }
        else // Skipping unsegregated runs only if the skipUnsegregatedData flag is true
            printf("Run %i did not show segregation; skipping the run.\n", runIterator);
            
        FreeAllocatedMemory(&read_data);
        runIterator++;
    } while(runIterator <= numberOfCalculationRuns);

    printf("Debugging: Number of snapshots that contributed to distribution: %i\n", counter);

    // Normalizing to get a probability distribution:
    CalculateDistribution(&dist_data, distribution, distributionLength);

    // Printing the probability distribution; if free energy is not being printed
    if(!calculateFreeEnergy)
    {
        // Getting free energy label:
        char* freeEnergyLabel; // A label indicating which reaction coordinate is being used
        SetFreeEnergyLabel(&freeEnergyLabel, acceptedReactionCoords[chosenReactionCoord], region1, region2);
        PrintProbabilityDistribution(&dist_data, distribution, freeEnergyLabel);
        free(freeEnergyLabel);
    }

    // Writing the successful exit message:
    char* regionLabel;
    if(DO_NOT_USE_REGIONS)
        regionLabel = "polymers";
    else
        regionLabel = "regions";
    printf("Computed probability distribution for the reaction coordinate '%s' between %s %i and %i using %li snapshots of the system.\n", acceptedReactionCoords[chosenReactionCoord], regionLabel, region1, region2, dist_data.numberOfDataPoints);

    FreeMonomerDistribution(&dist_data);
    
}

// Calculating and printing the free energy for each bin of the reaction coordinate using the probaility distribution
void ComputeFreeEnergy(int distributionLength, double distribution[distributionLength], char* freeEnergyLabel)
{
    double minimumProbability = 0.000000000001; // The minimum probability below which any probability will be considered as 0

    double freeEnergy[distributionLength];
    for(int i = 0; i < distributionLength; i++)
    {
        double probability = distribution[i] * dist_data.binWidth; // calculating probability from the probability density
        // checking if the probability is too low; close to 0:
        if(probability < minimumProbability)
            freeEnergy[i] = INFINITY; // will this work for plotting?
        else
            freeEnergy[i] = - 1 * log(probability); // F(x) = - kT ln (p(x)) upto an additive constant
    }

    // Printing the free energy with the reaction coordinate:
    char* writeFilePath;
    SetFreeEnergyFilePath(&writeFilePath, runIndex, freeEnergyLabel);
    FILE* writeFile = fopen(writeFilePath, "w");
    if(writeFile == NULL)
    {
        printf("ERROR: The free energy file %s could not be opened!\n", writeFilePath);
        exit(1);
    }
    // Writing header:
    fprintf(writeFile, "Number of simulation snapshots considered: %li\n", dist_data.numberOfDataPoints);
    fprintf(writeFile, "Reaction Coordinate, Free Energy (kT)\n"); 
    // Writing free energy:
    for(int i = 0; i < distributionLength; i++)
    {
        double reactionCoordinate = lowerBoundary + i * dist_data.binWidth;
        fprintf(writeFile, "%lf, %lf\n", reactionCoordinate, freeEnergy[i]);
    }
    fclose(writeFile);
    free(writeFilePath);

    printf("Debugging: Printed the free energies for the reaction coordinate using %li snapshots of the system.\n", dist_data.numberOfDataPoints);
}

int main(int argc, char** argv)
{
    SetConstants(argc, argv);

    // Initializing and reading the segregation times, if required:
    if(skipUnsegregatedData)
        InitializeSegregationTimes(&segData, numberOfMonomers, architecture, numberOfRuns);
    // PrintSegregationTimes(&segData); debugging

    // Initializing region:
    // InitializeRegionManually(&regData, numberOfPolymers, numberOfMonomers, true, false);
    InitializeRegion(&regData, numberOfPolymers, numberOfMonomers);

    // Computing free energy for multiple reaction coordinates: CoM distance between different regions
    // The region pairs have been defined in regionPairs
    for(int n = 0; n < numberOfPairs; n++)
    {
        char* freeEnergyLabel;
        SetFreeEnergyLabel(&freeEnergyLabel, acceptedReactionCoords[chosenReactionCoord], regionPairs[n][0], regionPairs[n][1]);
        // Printing the probability distribution for the reaction coordinate: z CoM distance
        double distribution[numberOfBins];
        ComputeProbabilityDistribution(numberOfBins, distribution, regionPairs[n][0], regionPairs[n][1]);
        if(calculateFreeEnergy)
            ComputeFreeEnergy(numberOfBins, distribution, freeEnergyLabel);
        free(freeEnergyLabel);
    }

    FreeRegionPairs();
    FreeSegregationTimes(&segData);
    FreeRegionMemory(&regData);
    return 0;
}