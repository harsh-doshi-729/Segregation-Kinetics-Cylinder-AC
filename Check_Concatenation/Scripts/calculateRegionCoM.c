// This script calculates the CoM of each region as defined by the regions_config.h for all timesteps 
// of a LAMMPS simulation position dump

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>

// importing paths file:
#include "../../../../../Scripts/System_File_Paths/system_file_paths.h"
#include GENERAL_FILE_METHODS
#include LAMMPS_POSITION_FILE
#include REGION

// Global variables:
int numberOfMonomers; // The number of monomers in a single polymer
int numberOfPolymers = 2; // The number of polymers in the system
int totalMonomers; // The total number of monomers in the system
char* architecture;

// Run information:
int runIndex;

// File information:
char* directory; // The path to the directory where all the concatenation dump files are stored
simulation_data read_data; // The struct for the LAMMPS Position dump file

// region data:
region_data regionData; // The struct to store and read the information about regions in a polymer

// Accepts the arguments from the command line and sets them to global variables
void SetConstants(int argc, char** argv)
{
    int numberOfMandatoryArguments = 3;
    if(argc < numberOfMandatoryArguments + 1)
    {
        printf("Not enough arguments passed! Please pass the number of monomers (integer), architecture(string), and runNumber (integer) as arguments from the command line\n");
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
        SetDirectoryPath(&directory, CHECK_CONCATENATION, numberOfMonomers);
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

// Sets the write file path for the CoM file of a region
// Expects the region number starting from 1
void SetWriteFilePath(char** filePathPointer, int regionNumber)
{
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    char* fileName;
    int bytes = asprintf(&fileName, "com%i.dat", regionNumber);
    if(bytes == -1)
    {
        printf("Memory could not be allocated for the the filename of the write file! Terminating.\n");
        exit(1);
    }
    Append(filePathPointer, directory, folder, fileName);
    free(folder);
    free(fileName);
}

// Calculates and writes the CoM for the current timestep
void CalculateCurrentCoM(FILE** writeFilePointers)
{
    // Array for CoMs:
    double CoM_Array[regionData.totalNumberOfRegions][3];
    // Initializing:
    for(int i = 0; i < regionData.totalNumberOfRegions; i++)
    {
        for(int j = 0; j < 3; j++) // iterating over the coordinates
            CoM_Array[i][j] = 0;
    }

    // Reading and Computing CoM:
    for(int n = 0; n < read_data.numberOfMonomers; n++)
    {
        int regionIndex = GetRegionID(regionData, n+1);
        for(int j = 0; j < 3; j++)
            CoM_Array[regionIndex][j] += read_data.bead_positions[n][j];
    }
    // Normalizing and Writing:
    for(int i = 0; i < regionData.totalNumberOfRegions; i++)
    {
        // Normalizing with number of monomers in the region:
        int regionCount = GetRegionCount(regionData, i);
        for(int j = 0; j < 3; j++)
            CoM_Array[i][j] = CoM_Array[i][j] / regionCount;
        
        // Writing to separate files:
        fprintf(writeFilePointers[i], "%i %lf %lf %lf\n", read_data.timeStep, CoM_Array[i][0], CoM_Array[i][1], CoM_Array[i][2]);
    }
    
}

// Calculates and writes the CoM of each timestep one at a time
void CalculateCoM(void)
{
    // Initializing Reading struct:
    char* readFilePath;
    SetReadFilePath(&readFilePath);
    InitializeSimulationData(&read_data, totalMonomers, readFilePath);
    free(readFilePath);
    // Initializing Region struct:
    InitializeRegion(&regionData, numberOfPolymers, numberOfMonomers);

    // Opening write files:
    FILE* writeFilePointers[regionData.totalNumberOfRegions];
    for(int i = 0; i < regionData.totalNumberOfRegions; i++)
    {
        char* writeFilePath;
        SetWriteFilePath(&writeFilePath, i+1);
        writeFilePointers[i] = fopen(writeFilePath, "w");
        if(writeFilePointers[i] == NULL)
        {
            printf("The write file '%s' could not be opened! Terminating.\n", writeFilePath);
            exit(1);
        }
        free(writeFilePath);
        // Writing file headers:
        fprintf(writeFilePointers[i], "Timestep, CoM_x, CoM_y, CoM_z\n");
    }

    // Reading, Calculating CoM, and Writing:
    while(ReadPositions(&read_data))
    {
        CalculateCurrentCoM(writeFilePointers);
    }

    // Closing files:
    for(int i = 0; i < regionData.totalNumberOfRegions; i++)
    {
        fclose(writeFilePointers[i]);
    }

    // Freeing memory:
    FreeAllocatedMemory(&read_data);
    FreeRegionMemory(&regionData);
}

int main(int argc, char** argv)
{
    SetConstants(argc, argv);
    CalculateCoM();
    free(directory);
    return 0;
}