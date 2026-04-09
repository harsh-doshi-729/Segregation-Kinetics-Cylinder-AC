// This script calculates the monomer density in sections of a system of polymers by using the header files LAMMPS_Position_File.h and MonomerDistribution.h

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>

#include "../../Global_Scripts/System_File_Paths/system_file_paths.h"
#include GENERAL_FILE_METHODS
#include LAMMPS_POSITION_FILE
#include MONOMER_DISTRIBUTION
#include SECTION
#include CROSS_LINKS_DATABASE

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
int equilibrationSteps; // The number of iterations/timesteps in the simulation when equilibration occurs


// Reading positions:
char* directory; // The path to the directory where the simulation files are stored for the desired run
char* initializationProcedure; // The name of the initialization procedure used to generate the initial state; this is used to read from the correct folder
simulation_data read_data; // The struct to read the data

// Calculating monomer distribution:
double binWidth = 0.4;
// int numberOfSections = 10; // The total number of sections in the cylinder is divided in
monomer_distribution_data* dist_array; // An array of monomer_dsitribution_data structs for each of the different sections in the system

// The architecture cross-links object for reading confinement diameters
Architecture archDiameter;

// Reads the arguments from the command line and sets the corresponding constants:
void SetConstants(int argc, char** argv)
{
    int numberOfMandatoryArguments = 5;
    if(argc < numberOfMandatoryArguments + 1)
    {
        printf("Not enough arguments passed! Please pass the number of monomers (integer), architecture(string), run number (integer), the bin width (float), and the initialization procedure (string) as arguments from the command line in the following format:\n");
        printf("./sectional.out <numberOfMonomers> <architecture> <runIndex> <binWidth> <initializationProcedure>\n");
        printf("Optionally an argument for the finite axis length (float) can be passed after the %i arguments.\n", numberOfMandatoryArguments);
        exit(1);
    }
    else
    {
        numberOfMonomers = atoi(argv[1]); // the first argument after the program name
        architecture = argv[2]; // the second argument after the program name
        runIndex = atoi(argv[3]); // The simulation data files will be in the run1 folder by default
        binWidth = atof(argv[4]);
        initializationProcedure = argv[5]; // the fifth argument after the program name

        // Verifying the arguments:
        if(numberOfMonomers == 0)
        {
            printf("The argument passed for number of monomers %s could not be converted to an integer! Terminating.\n", argv[1]);
            exit(1);
        }
        if(binWidth <= 0)
        {
            printf("The bin width must be positive! Terminating.\n");
            exit(1);
        }
        

        totalMonomers = numberOfMonomers * numberOfPolymers;
        containsPreEquilibriumData = false;
        equilibrationSteps = 2 * pow(10, 8);

        char* folderPrefixPath;
        SetAbsolutePath(&folderPrefixPath, CREATE_INITIAL_STATES);
        SetDirectoryPath(&directory, folderPrefixPath, numberOfMonomers, initializationProcedure);

        // Default value of axis length
        isLengthInfinite = true;
        if(numberOfMonomers == 200)
        {
            axisLength = 500; // Length of long cylinder in case of 200 monomers
        }
        else if(numberOfMonomers == 500)
        {
            axisLength = 1000; // Length of long cylinder in case of 500 monomers
        }
        // Initializing architecture object
        InitializeArcAndReadDiameter(&archDiameter, numberOfMonomers, architecture, initializationProcedure);

        // Checking for optional argument
        if(argc > numberOfMandatoryArguments + 1)
        {
            isLengthInfinite = false;
            axisLength = atof(argv[numberOfMandatoryArguments+1]);
            if(axisLength <= 0)
            {
                printf("The axis length must be positive! Terminating.\n");
                exit(1);
            }
        }
    }
}

// Allocates memory for the monomer distribution structs array
void AllocateDistributionStructArray(void)
{
    dist_array = malloc(numberOfSections * sizeof(*dist_array));
    assert(dist_array);
}

// Sets the file path for the position dump file to be read from
void SetReadFilePath(char** filePathPointer)
{
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    char* fileName = "distribution_positions.dump";
    Append(filePathPointer, directory, folder, fileName);
    free(folder);
}

// Sets the path to the file where the monomer distribution for a particular region is to be written
void SetWriteFilePath(char** filePathPointer, int sectionID)
{
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    char* fileName; // with the preceding folder too
    int bytes = asprintf(&fileName, "monomer_distribution/monomer_distribution_sect%i.csv", sectionID);
    if(bytes == -1)
    {
        printf("Memory could not be allocated for the write file name! Terminating.\n");
        exit(1);
    }
    Append(filePathPointer, directory, folder, fileName);
    free(folder);
    free(fileName);
}

// Returns the total number of data points used to generate the distributions of all sections
int GetTotalNumberOfDataPoints(void)
{
    int numberOfDataPoints = 0;
    for(int n = 0; n < numberOfSections; n++)
    {
        numberOfDataPoints += dist_array[n].numberOfDataPoints;
    }
    return numberOfDataPoints;
}

// Reads the data file, computes the distribution for each section. and writes them to a file
void ComputeEntireDistribution(void)
{
    // Reading data file:
    char* readFilePath;
    SetReadFilePath(&readFilePath);
    // Note: The read_data struct has already been allocated static memory while declaration
    InitializeSimulationData(&read_data, totalMonomers, readFilePath);

    // Initializing monomer distribution structs:
    AllocateDistributionStructArray();
    for(int n = 0; n < numberOfSections; n++)
    {
        InitializeDistributionDataWithAxisLengthBinWidth(&(dist_array[n]), totalMonomers, axisLength, isBoxCentredAtOrigin, binWidth);
    }

    double radius = archDiameter.confinementDiameter / 2;

    // Reading File:
    while(ReadPositions(&read_data))
    {
        // Skipping pre equilibrium data:
        if(containsPreEquilibriumData && read_data.timeStep < equilibrationSteps)
            continue;
        
        // Adding each monomer to its region distribution:
        for(int i = 0; i < totalMonomers; i++)
        {
            int sectionID = GetSectionID(numberOfMonomers, radius, numberOfSections, read_data.bead_positions[i]);
            IncrementMonomerCounts(&(dist_array[sectionID - 1]), read_data.bead_positions[i][2]);
        }
    }
    CloseReadFile(&read_data);
    FreeAllocatedMemory(&read_data);
    free(readFilePath);
    FreeArchitecture(&archDiameter); // Done using the confinement diameter information
    
    // Computing and writing the distribution to file:
    int normalizationFactor = GetTotalNumberOfDataPoints();
    for(int n = 0; n < numberOfSections; n++)
    {
        char* writeFilePath;
        SetWriteFilePath(&writeFilePath, n+1);
        FILE* writeFilePointer = fopen(writeFilePath, "w");
        if(writeFilePointer == NULL)
        {
            printf("The write file at %s could not be opened! Terminating.\n", writeFilePath);
            exit(1);
        }
        double distribution[dist_array[n].numberOfBins]; // The distribution array to be calculated; will be passed as is
        CalculateDistributionWithNormalizationFactor(&dist_array[n], distribution, dist_array[n].numberOfBins, normalizationFactor);
        PrintDistribution(&(dist_array[n]), distribution, isLengthInfinite, writeFilePointer);
        fclose(writeFilePointer);
        free(writeFilePath);
        FreeMonomerDistribution(&(dist_array[n]));
    }
    free(dist_array); // freeing the monomer distribution struct array
}

int main(int argc, char** argv)
{
    // Testing the GetSectionID method:
    // double position[3] = {0.6839, 1.8603, -11.8};
    // printf("Section ID: %i\n", GetSectionID(position));

    printf("Number of Sections: %i\n", numberOfSections); // Printing the number of sections as set in Section.h
    SetConstants(argc, argv);
    ComputeEntireDistribution();
    printf("The monomer distributions for different sections have been written for %s Run %i.\n", architecture, runIndex);
}
