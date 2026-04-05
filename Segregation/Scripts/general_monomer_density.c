// This script calculates the monomer density of a system of polymers by using the header files LAMMPS_Position_File.h and MonomerDistribution.h

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
#include MONOMER_DISTRIBUTION
#include REGION

// Global variable:
int numberOfMonomers; // The number of monomers in a single polymer
int numberOfPolymers = 2; // The number of polymers in the system
int totalMonomers; // The total number of monomers in the system
char* architecture;
double axisLength; // The length of the axis along which the distribution needs to be calculated
bool isLengthInfinite; // A flag to indicate whether the length of the cylinder is very long as compared to the extent of the polymers
bool isBoxCentredAtOrigin = true;

// Run information:
int runIndex;
bool containsPreEquilibriumData; // A flag to indicate whether the simulation dump file contains data taken before equilibration
bool imposeUpperBound; // A flag to indicate whether the monomer density should only be calculated till some upper bound in timesteps
int equilibrationSteps; // The number of iterations/timesteps in the simulation when equilibration occurs
int timeUpperBound;
bool isSystemUncut; // A flag to indicate whether the system is "uncut" i.e. a control for the a system of polymers with their loops cut. 
// For systems that have nothing to do with LoopCutting at all, set this above flag to false.

// Reading positions:
char* initializationProcedure; // The name of the procedure used to generate the mixed state; for example: "fene_recenter", "fene_glued", etc.
char* destinationLabel; // A string to store a label for the destination. Currently accepted labels: "new_segregation" and "Create_Initial_States"
char* acceptedDestinations[2] = {"new_segregation", "Create_Initial_States"};
char* directory; // The path to the directory where the simulation files are stored for the desired run
simulation_data read_data; // The struct to read the data

// region data:
region_data regionData; // The struct to store and read the information about regions in a polymer

// Calculating monomer distribution:
double binWidth = 0.2;
monomer_distribution_data* dist_array; // An array of monomer_dsitribution_data structs for each of the different regions in the system

// Calculation information:
bool useSingleSnapshot = false; // A flag indicating whether only a single snapshot should be used to calculate the distribution
// If true, the last snapshot will be used to calculate the distribution

// Reads the arguments from the command line and sets the corresponding constants:
void SetConstants(int argc, char** argv)
{
    int numberOfMandatoryArguments = 5;
    if(argc < numberOfMandatoryArguments + 1) // +1 since the script name is passed as the first argument always
    {
        printf("Not enough arguments passed! Please pass the number of monomers (integer), architecture(string), runNumber(integer), destination label (string), and the name of the initialization procedure (string) as arguments from the command line\n");
        printf("Accepted destination labels: ");
        PrintArray(2, acceptedDestinations);
        printf("Optional arguments for a boolean argument (0 or 1) to indicate if a single snapshot should be used to compute the distribution and the length of the cylinder (double) can be passed after the above arguments.\n");
        exit(1);
    }
    else
    {
        numberOfMonomers = atoi(argv[1]); // the first argument after the program name
        architecture = argv[2]; // the second argument after the program name
        runIndex = atoi(argv[3]); // the third argument after the program name
        destinationLabel = argv[4];
        initializationProcedure = argv[5];
        isLengthInfinite = true; // Default value in case the axis length is not passed
        isSystemUncut = false; // Default value for the system being uncut; can be changed by passing the appropriate argument from the command line
        // Verifying the arguments:
        if(numberOfMonomers == 0)
        {
            printf("The argument passed for number of monomers %s could not be converted to an integer! Terminating.\n", argv[1]);
            exit(1);
        }
        if(runIndex == 0)
        {
            printf("The argument passed for the run index %s could not be converted to an integer! Terminating.\n", argv[3]);
            exit(1);
        }

        totalMonomers = numberOfMonomers * numberOfPolymers;
        containsPreEquilibriumData = false;
        if(numberOfMonomers == 200)
        {
            equilibrationSteps = 4 * pow(10, 6); // single relaxation time
            axisLength = 500; // Length of long cylinder in case of 200 monomers
        }
        else if(numberOfMonomers == 500)
        {
            equilibrationSteps = 2.5 * pow(10, 7); // single relaxation time
            axisLength = 1000; // Length of long cylinder in case of 500 monomers
        }
        imposeUpperBound = false;
        timeUpperBound = equilibrationSteps; // reusing equilibration steps as the upper bound for timesteps
        if(imposeUpperBound)
        {
            printf("Note: The monomer density will only be calculated upto the timestep %i! To remove this bound, edit the 'imposeUpperBound' flag in the C script.\n", timeUpperBound);
        }

        char* folderPrefix;
        if(strcasecmp(destinationLabel, acceptedDestinations[0]) == 0) // new_segregation chosen
            folderPrefix = SEGREGATION;
        else if(strcasecmp(destinationLabel, acceptedDestinations[1]) == 0) // Create_Initial_States chosen
            folderPrefix = CREATE_INITIAL_STATES;
        else
        {
            printf("The destination label passed %s is not recognized! Please pass one form the following options: ", destinationLabel);
            PrintArray(2, acceptedDestinations);
            exit(1);
        }
        char* prefixAbsDirectory; // absolute path to the folder prefix directory
        SetAbsolutePath(&prefixAbsDirectory, folderPrefix);
        SetDirectoryPath(&directory, prefixAbsDirectory, numberOfMonomers, initializationProcedure);
        free(prefixAbsDirectory);

        // Optional arguments:
        if(argc > numberOfMandatoryArguments + 1)
        {
            // The first optional argument is assumed to be the boolean flag
            if(atoi(argv[numberOfMandatoryArguments+1]) == 1) // If one was passed
            {
                useSingleSnapshot = true;
                printf("Computing the distribution from a single snapshot.\n");
            }
            else
                printf("Ignoring the optional argument %s since it seems meaningless.\n", argv[numberOfMandatoryArguments+1]);

            if(argc == numberOfMandatoryArguments + 3) // Two optional arguments passed
            {
                axisLength = atof(argv[numberOfMandatoryArguments+2]);
                isLengthInfinite = false;
            }
        }

        // Debugging:
        if(isLengthInfinite)
        {
            printf("No axis length passed. Assuming cylinder to be infinite in length.\n");
        }
        else
        {
            printf("The axis length is set to %lf.\n", axisLength);
        }

    }
}

// Allocates memory for the monomer distribution structs array
void AllocateDistributionStructArray(void)
{
    dist_array = malloc(regionData.totalNumberOfRegions * sizeof(*dist_array));
    assert(dist_array);
}

// Sets the file path for the position dump file to be read from
void SetReadFilePath(char** filePathPointer)
{
    char* folder;
    SetUncutFolderName(&folder, architecture, runIndex, isSystemUncut);
    char* fileName;
    if(strcasecmp(destinationLabel, acceptedDestinations[0]) == 0) // new_segregation chosen
        fileName = "visual.dump";
    else if(strcasecmp(destinationLabel, acceptedDestinations[1]) == 0) // Create_Initial_States chosen
        fileName = "distribution_positions.dump";

    Append(filePathPointer, directory, folder, fileName);
    free(folder);
}

// Sets the path to the file where the monomer distribution for a particular region is to be written
// Assuming the regions are indexed from 1
void SetWriteFilePath(char** filePathPointer, int regionID)
{
    char* folder;
    SetUncutFolderName(&folder, architecture, runIndex, isSystemUncut);
    char* fileNamePrefix;
    if(strcasecmp(destinationLabel, acceptedDestinations[0]) == 0) // new_segregation chosen
        fileNamePrefix = "monomer_distribution_reg";
    else if(strcasecmp(destinationLabel, acceptedDestinations[1]) == 0) // Create_Initial_States chosen
        fileNamePrefix = "monomer_distribution/monomer_distribution_reg";
    
    char* singleSnapshotSuffix = "";
    if(useSingleSnapshot)
        singleSnapshotSuffix = "_single"; // A suffix used when only a single snapshot is used to calculate the distribution
    char* fileName;
    int bytes = asprintf(&fileName, "%s%i%s.csv", fileNamePrefix, regionID, singleSnapshotSuffix);
    if(bytes == -1)
    {
        printf("Memory could not be allocated for the write file name! Terminating.\n");
        exit(1);
    }

    Append(filePathPointer, directory, folder, fileName);
    free(folder);
    free(fileName);
}

// Allocates memory to and initializes various structs: simulation_data, region_data, monomer_dist_data
void AllocateAndInitializeStructs(char* readFilePath)
{
    // Initializing reading machinery:
    // Note: The read_data struct has already been allocated static memory during declaration
    InitializeSimulationData(&read_data, totalMonomers, readFilePath);

    // Initializing and allocating the region_data struct:
    InitializeRegion(&regionData, numberOfPolymers, numberOfMonomers);

    // Initializing monomer distribution structs:
    AllocateDistributionStructArray();
    for(int n = 0; n < regionData.totalNumberOfRegions; n++)
    {
        InitializeDistributionDataWithAxisLengthBinWidth(&(dist_array[n]), GetNumberOfMonomersInRegion(regionData, n), axisLength, isBoxCentredAtOrigin, binWidth);
    }
}

// Computes and prints the distribution after the position file is read
void ComputeAndPrintDistribution(void)
{
    // Computing and writing the distribution to file:
    int totalDataPoints = 0;
    for(int n = 0; n < regionData.totalNumberOfRegions; n++)
        totalDataPoints += dist_array[n].numberOfDataPoints;

    for(int n = 0; n < regionData.totalNumberOfRegions; n++)
    {
        char* writeFilePath;
        SetWriteFilePath(&writeFilePath, n+1);
        FILE* writeFilePointer = fopen(writeFilePath, "w");
        if(writeFilePointer == NULL)
        {
            printf("The file %s could not be opened for writing! Terminating.\n", writeFilePath);
            exit(1);
        }
        
        // calculating distribution:
        int arraySize = dist_array[n].numberOfBins;
        double distribution[arraySize];
        // Calculating the distribution with a global normalizing factor of the total number of data points for all monomers:
        // Adding up the probability distributions of all regions will give the normalized total probability distribution
        CalculateDistributionWithNormalizationFactor(&(dist_array[n]), distribution, arraySize, totalDataPoints);
        PrintDistribution(&(dist_array[n]), distribution, isLengthInfinite, writeFilePointer);
        fclose(writeFilePointer);
        free(writeFilePath);
    }
}

// Computes the distribution of all regions over the entire simulation
void ComputeEntireDistribution(void)
{
    char* readFilePath;
    SetReadFilePath(&readFilePath);
    AllocateAndInitializeStructs(readFilePath);
    free(readFilePath); // The struct contains a copy of the file path; this one is redundant now

    // Reading File:
    while(ReadPositions(&read_data))
    {
        // Skipping pre equilibrium data:
        if(containsPreEquilibriumData && read_data.timeStep < equilibrationSteps)
            continue;
        
        // imposing upper bound for timesteps
        if(imposeUpperBound && read_data.timeStep > timeUpperBound)
        {
            printf("Timestep Upper Bound Reached.\n");
            break;
        }

        // Adding each monomer to its region distribution:
        for(int i = 0; i < totalMonomers; i++)
        {
            int regionID = GetRegionID(regionData, i+1); // regionID is 0 based
            IncrementMonomerCounts(&(dist_array[regionID]), read_data.bead_positions[i][2]);
        }
    }
    CloseReadFile(&read_data);
    FreeAllocatedMemory(&read_data);
    
    ComputeAndPrintDistribution();
    // freeing struct memory:
    for(int n = 0; n < regionData.numberOfRegions; n++)
        FreeMonomerDistribution(&(dist_array[n]));
    free(dist_array); // freeing the monomer distribution struct array
    FreeRegionMemory(&regionData);
}

// Computes the monomer distribution based on a single snapshot: the last snapshot of the read file
void ComputeSingleSnapshotDistribution(void)
{
    char* readFilePath;
    SetReadFilePath(&readFilePath);
    AllocateAndInitializeStructs(readFilePath);
    free(readFilePath); // The struct contains a copy of the file path; this one is redundant now

    // Reading the file:
    while(ReadPositions(&read_data)); // iterating till EOF is reached
    // The last snapshot should be stored in the read_data struct
    // Adding each monomer to its region distribution:
    for(int i = 0; i < totalMonomers; i++)
    {
        int regionID = GetRegionID(regionData, i+1); // regionID is 0 based
        IncrementMonomerCounts(&(dist_array[regionID]), read_data.bead_positions[i][2]);
    }
    CloseReadFile(&read_data);
    FreeAllocatedMemory(&read_data);
    
    ComputeAndPrintDistribution();
    // freeing struct memory:
    for(int n = 0; n < regionData.numberOfRegions; n++)
        FreeMonomerDistribution(&(dist_array[n]));
    free(dist_array); // freeing the monomer distribution struct array
    FreeRegionMemory(&regionData);
}

int main(int argc, char** argv)
{
    SetConstants(argc, argv);
    if(!useSingleSnapshot)
        ComputeEntireDistribution();
    else
        ComputeSingleSnapshotDistribution();
}