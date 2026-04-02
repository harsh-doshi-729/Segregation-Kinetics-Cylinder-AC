// This script is used to calculate the Centres of Mass of various regions of polymers at each timestep
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>

// importing paths file:
#include "../../../../Scripts/System_File_Paths/system_file_paths.h"
#include  GENERAL_FILE_METHODS
#include LAMMPS_POSITION_FILE
#include REGION

// Global variables:

// System information:
int numberOfMonomers; // The number of monomers in a single polymer
int numberOfPolymers = 2; // The number of polymers in the system
int totalMonomers; // The total number of monomers in the system
char* architecture; // The name of the architecture

// Simulation information:
int runIndex;

// Reading file tools:
char* directory; // The path to the folder where the simulation data is stored
simulation_data read_data; // The struct to read a LAMMPS Position fump file

// Region information:
region_data regionData; // The struct that stores region information

// Options for the directory:
    int arrayLength = 2;
    char* acceptedDirectoryLabels[] = {"new_segregation", "Create_Initial_States"};
    char* acceptedDirectories[] = {NEW_SEGREGATION, CREATE_INITIAL_STATES}; // The corresponding directory paths
    char* directoryLabel; // The label for the desired directory
    int directoryIndex; // The index in the directory arrays for the desired directory

void SetConstants(int argc, char** argv)
{
    directoryIndex = 0;
    directoryLabel = acceptedDirectoryLabels[directoryIndex]; // optional argument; default value: new_segregation
    int numberOfMandatoryArguments = 3;


    if(argc < numberOfMandatoryArguments + 1)
    {
        printf("Not enough arguments passed! Please pass the number of monomers, the name of the architecture, and the run index.\n");
        printf("Format: ./a.out <numberOfMonomers> <architecture> <runIndex>\n");
        printf("An optional argument for the directory can be passed. Accepted options: ");
        PrintArray(arrayLength, acceptedDirectoryLabels);
        exit(1);
    }
    else
    {
        numberOfMonomers = atoi(argv[1]);
        // TODO: Make common functions for verification of the passed arguments
        if(numberOfMonomers == 0)
        {
            printf("The argument passed for number of monomers %s could not be converted to an integer! Terminating.\n", argv[1]);
            exit(1);
        }

        architecture = argv[2];

        runIndex = atoi(argv[3]);
        if(runIndex == 0)
        {
            printf("The argument passed for run index %s could not be converted to an integer! Terminating.\n", argv[3]);
            exit(1);
        }

        totalMonomers = numberOfPolymers * numberOfMonomers;

        if(argc > numberOfMandatoryArguments + 1) // optional argument passed
        {
            directoryLabel = argv[numberOfMandatoryArguments+1];
        }

        // Setting directory:
        bool isArgumentValid = false; // A flag to indicate whether the destination label argument is valid
        for(int i = 0; i < arrayLength; i++)
        {
            // Comparing the argument with accepted directories:
            if(strcasecmp(directoryLabel, acceptedDirectoryLabels[i]) == 0)
            {
                isArgumentValid = true; // since the arguments matched one of the acceptable options
                SetDirectoryPath(&directory, acceptedDirectories[i], numberOfMonomers);
                directoryIndex = i;
                break;
            }
        }
        if(!isArgumentValid)
        {
            printf("The argument %s was not recognized! Please enter one of the accepted directory labels: ", directoryLabel);
            PrintArray(arrayLength, acceptedDirectoryLabels);
            exit(1);
        }
    }
}

// Sets the file path to the position dump file to be read from
void SetReadFilePath(char** filePathPointer)
{
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    char* fileNames[] = {"visual.dump", "distribution_positions.dump"};
    char* fileName = fileNames[directoryIndex];
    Append(filePathPointer, directory, folder, fileName);
    free(folder);
}

// Sets the file path where the CoM time series for a particular region will be written;
// regionNumber >= 1
void SetWriteFilePath(char** filePathPointer, int regionNumber)
{
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    char* fileName;
    // Setting filename
    int bytes = asprintf(&fileName, "com_reg%i.dat", regionNumber);
    if(bytes == -1)
    {
        printf("Memory could not be allocated for the write file path! Terminating.\n");
        exit(1);
    }

    Append(filePathPointer, directory, folder, fileName);
    // Freeing folder and filename:
    free(folder);
    free(fileName);
}

// Opens the write files for each region's CoM
void OpenWriteFiles(FILE** writeFilePointers)
{
    for(int n = 0; n < regionData.totalNumberOfRegions; n++)
    {
        char* writeFilePath;
        SetWriteFilePath(&writeFilePath, n + 1);
        writeFilePointers[n] = fopen(writeFilePath, "w");
        if(writeFilePointers[n] == NULL)
        {
            printf("The write file %s could not be opened! Terminating.\n", writeFilePath);
            exit(1);
        }
        free(writeFilePath);
    }
}

// Prints the headers for the region CoM time series files
void PrintFileHeaders(FILE** writeFilePointers)
{
    for(int n = 0; n < regionData.totalNumberOfRegions; n++)
    {
        fprintf(writeFilePointers[n], "Time Step, x CoM, y CoM, z CoM\n");
    }
}

// Initializes the values of a 2D double array to 0s (actually a look up table)
void Initialize2DArray(double*** arrayPointer, size_t x, size_t y)
{
    *arrayPointer = (double**) malloc(x * sizeof(double*));
    for(int i = 0; i < x; i++)
    {
        (*arrayPointer)[i] = (double*) calloc(y, sizeof(double));
    }
}

// Frees the memory allocated for the 2D array
void Free2DArray(double** array, size_t x, size_t y)
{
    for(int i = 0; i < x; i++)
        free(array[i]);
    free(array);
}

// Initializes the values of an int array to 0s
void InitializeArray(int* array, size_t arrayLength)
{
    for(int i = 0; i < arrayLength; i++)
        array[i] = 0;
}

// Calculates and writes the CoMs of all regions for a single timestep
void PrintSingleTimestepCoM(FILE** writeFilePointers)
{
    // Defining a few arrays to individually calculate CoM for each region:
    double** CoM_Array; // 2D array; 1st dim = regionID; 2nd dim = spatial coordinate
    int counter_Array[regionData.totalNumberOfRegions]; // An array of counters to count the number of monomers in each region

    Initialize2DArray(&CoM_Array, regionData.totalNumberOfRegions, 3);
    InitializeArray(counter_Array, regionData.totalNumberOfRegions);

    // Assuming the positions for the current timestep have already been read
    // Assuming the region data struct has been initialized
    // Assuming all monomers have the same mass

    //debugging:
    // int regionIDs[GetRegionCount(&regionData, 0)];
    // int counter = 0;
    // debugging: printing the monomer index and region ID headers:
    // printf("Monomer Index, RegionID\n");
    for(int i = 0; i < read_data.numberOfMonomers; i++)
    {
        int regionID = GetRegionID(regionData, i + 1);
        // debugging: printing the monomer index and region ID for each monomer:
        // printf("%i, %i\n", i+1, regionID);
        // debugging: checking region IDs of all monomers that should be present in the first region
        // if(i+1 <= numberOfMonomers && (i+1 >= regionData.regions[0].lower || i+1 <= regionData.regions[0].upper))
            // regionIDs[counter++] = regionID;

        // Adding the coordinates of the monomer to the relevant regions's CoM:
        for(int j = 0; j < 3; j++) // iterating over spatial dimensions
            CoM_Array[regionID][j] += read_data.bead_positions[i][j];
        counter_Array[regionID]++; // incrementing the region counter
    }

    // Averaging each region's CoM:
    for(int n = 0; n < regionData.totalNumberOfRegions; n++)
    {
        // Verifying that the number of monomers encountered matches the expected region count:
        int regionCount = GetRegionCount(regionData, n);
        if(counter_Array[n] != regionCount)
        {
            printf("The number of monomers counted %i does not match the expected region count %i for region %i! Terminating.\n", counter_Array[n], regionCount, n + 1);
            exit(1);
        }
        // else: verification successful:
        for(int j = 0; j < 3; j++) // iterating over spatial dimensions
            CoM_Array[n][j] = CoM_Array[n][j] / regionCount;
    }

    // Writing to the region CoM Time series files:
    for(int n = 0; n < regionData.totalNumberOfRegions; n++)
    {
        fprintf(writeFilePointers[n], "%i, %.6lf, %.6lf, %.6lf\n", read_data.timeStep, CoM_Array[n][0], CoM_Array[n][1], CoM_Array[n][2]);
    }

    Free2DArray(CoM_Array, regionData.totalNumberOfRegions, 3);
}
// TODO: Also print the min and max position of a bead for a region; this will help visualise the extent of a region

// Closes the write files for each region
void CloseWriteFiles(FILE** writeFilePointers)
{
    for(int n = 0; n < regionData.totalNumberOfRegions; n++)
    {
        fclose(writeFilePointers[n]);
    }
}

// Calculates and prints the CoM for all regions for all the timesteps
void CalculateAndPrintCoM(void)
{
    // Initializing the read data struct
    char* readFilePath;
    SetReadFilePath(&readFilePath);
    InitializeSimulationData(&read_data, totalMonomers, readFilePath);

    // Initializing the region data struct
    InitializeRegion(&regionData, numberOfPolymers, numberOfMonomers);

    // Defining write file pointers:
    FILE* writeFilePointers[regionData.totalNumberOfRegions];
    OpenWriteFiles(writeFilePointers);

    // Writing to files:
    PrintFileHeaders(writeFilePointers);
    while(ReadPositions(&read_data))
    {
        PrintSingleTimestepCoM(writeFilePointers);
    }
    CloseWriteFiles(writeFilePointers);

    // freeing allocated memory:
    FreeAllocatedMemory(&read_data);
    FreeRegionMemory(&regionData); 
    free(directory);

    // Writing output message:
    printf("CoM for different regions calculated for %s Run %i.\n", architecture, runIndex);
}

int main(int argc, char** argv)
{
    SetConstants(argc, argv);
    CalculateAndPrintCoM();
    return 0;
}