// This script calculates the radial monomer density of a system of polymers by using the header files LAMMPS_Position_File.h and MonomerDistribution.h

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>

#include "../../../../Scripts/System_File_Paths/system_file_paths.h"
#include GENERAL_FILE_METHODS
#include LAMMPS_POSITION_FILE
#include MONOMER_DISTRIBUTION
#include CROSS_LINKS_DATABASE
#include REGION

// global variables:
int numberOfMonomers; // number of monomers in a single polymer
int numberOfPolymers = 2;
int totalMonomers;
char* architecture;
int runIndex;

bool containsPreEquilibriumData;
int equilibrationSteps;

// distribution information:
double binWidth;
int numberOfBins;
Architecture archDiameter; // Architecture object to read the confining diameter
double radius; 
monomer_distribution_data dist_data;

// File I/O information:
char* directory; // The path to the folder where the desired file is stored
simulation_data read_data;

// Directory labels:
int labelArrayLength = 2;
char* acceptedDirectoryLabels[] = {"new_segregation", "Create_Initial_States"};
bool forInitialization = false; // by default, the distribution is calculated in Create_Initial_States

// Reads the arguments from the command line and sets the corresponding constants:
void SetConstants(int argc, char** argv)
{
    int numberOfMandatoryArguments = 4;
    if(argc < numberOfMandatoryArguments + 1)
    {
        printf("Not enough arguments passed! Please pass the number of monomers (integer), architecture(string), the runNumber(integer), and the binWidth (float) as arguments from the command line in the following format:\n");
        printf("./radial.out <numberOfMonomers> <architecture> <runIndex> <binWidth>\n");
        printf("An optional argument for the directory label can be passed to calculate the distribution in that directory.\n");
        printf("Accepted directory labels: ");
        PrintArray(labelArrayLength, acceptedDirectoryLabels);
        exit(1);
    }
    else
    {
        numberOfMonomers = atoi(argv[1]); // the first argument after the program name
        architecture = argv[2]; // the second argument after the program name
        runIndex = atoi(argv[3]); // the third argument after the program name
        binWidth = atof(argv[4]);

        // Verifying inputs:
        if(numberOfMonomers == 0)
        {
            printf("The argument passed for the number of monomers %s could not be converted to an integer! Terminating.\n", argv[1]);
            exit(1);
        }
        if(runIndex == 0)
        {
            printf("The argument passed for the run index %s could not be converted to an integer! Terminating.\n", argv[3]);
            exit(1);
        }
        if(binWidth == 0)
        {
            printf("The argument passed for the bin width %s could not be converted to a float! Terminating.\n", argv[4]);
            exit(1);
        }

        InitializeArcAndReadDiameter(&archDiameter, numberOfMonomers, architecture);
        radius = archDiameter.confinementDiameter / 2;

        totalMonomers = numberOfMonomers * numberOfPolymers;
        // numberOfBins = 100;
        // binWidth = 0.1;

        containsPreEquilibriumData = false;
        // equilibrationSteps = 2 * pow(10, 6);

        // Checking for optional argument:
        char* folderPrefix = NEW_SEGREGATION; // default directory 
        if(argc > numberOfMandatoryArguments + 1)
        {
            char* label = argv[numberOfMandatoryArguments+1];
            int index = SearchStringArray(label, labelArrayLength, acceptedDirectoryLabels);
            if(index == -1)
            {
                printf("The label %s passed as the optional argument could not be recognized! Please pass on of the following accepted labels: ", label);
                PrintArray(labelArrayLength, acceptedDirectoryLabels);
                exit(1);
            }
            else if(index == 1) // new_segregation passed
            {
                folderPrefix = CREATE_INITIAL_STATES;
                forInitialization = true;
            }
        }
        // Setting the directory:
        SetDirectoryPath(&directory, folderPrefix, numberOfMonomers);
    }
}

// Sets the read file path to the position dump file to be read from:
void SetReadFilePath(char** filePathPointer)
{
    char* fileName = "visual.dump"; // default filename outside Create_Initial_States
    if(forInitialization)
        fileName = "distribution_positions.dump";
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    Append(filePathPointer, directory, folder, fileName);
}

// Sets the write file path where the radial distribution will be written
// Args:
// - filePathPointer: pointer to the string where the file path will be stored
// - useRegions: boolean variable to indicate whether the distribution is being calculated for different regions of the polymer
// - regionIndex: the index of the region for which the distribution is being calculated (0-indexed); only used if useRegions is true
void SetWriteFilePath(char** filePathPointer, bool useRegions, int regionIndex)
{
    char* fileNamePrefix = "radial_monomer_distribution";
    if(forInitialization)
        fileNamePrefix = "monomer_distribution/radial_monomer_distribution"; // filename with the preceeding folder
    char* regionLabel = "";
    if(useRegions)
    {
        int bytes = asprintf(&regionLabel, "_reg%i", regionIndex + 1); // regionIndex is 0-indexed; adding 1 for human readability
        if(bytes == -1)
        {
            printf("Memory could not be allocated for the file name region label!\n");
            exit(1);
        }
    }
    char* fileName;
    Append(&fileName, fileNamePrefix, regionLabel, ".csv");
    if(useRegions)
        free(regionLabel);
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    Append(filePathPointer, directory, folder, fileName);
    free(fileName);
    free(folder);
}

// Calculates and returns the radial distance of a position vector from the z-axis
double GetRadialDistance(double* position)
{
    double radialDistance = sqrt(pow(position[0], 2) + pow(position[1], 2));
    return radialDistance;
}

// Calculates the radial probability distribution by normalizing the monomer counts correctly
void CalculateRadialDistribution(monomer_distribution_data dist_data, size_t arraySize, double distribution[arraySize])
{
    // total normalization factor = 2 * PI * R_avg * binWidth * number of frames; R_avg = (R_out + R_in) / 2 = R_in + binWidth / 2
    double preNormalization = 2 * M_PI * binWidth * dist_data.numberOfDataPoints; // common factor for normalizing all bins
    for(int i = 0; i < arraySize; i++)
    {
        double averageRadius = i * binWidth + binWidth / 2;
        double normalization = preNormalization *  averageRadius;
        distribution[i] = dist_data.monomerCounts[i] / normalization;
    }
}

// Computes and prints the radial distribution to file:
void ComputeRadialDistribution(void)
{
    // Setting up reading machinery:
    char* readFilePath;
    SetReadFilePath(&readFilePath);
    InitializeSimulationData(&read_data, totalMonomers, readFilePath);

    // Setting up the monomer distribution struct:
    InitializeDistributionDataWithAxisLengthBinWidth(&dist_data, totalMonomers, radius, false, binWidth);
    // InitializeDistributionDataWithAxisLengthNumberOfBins(&dist_data, totalMonomers, radius, false, numberOfBins);
    // The radial distribution ranges from 0 to the radius

    // Reading data and adding to distribution:
    while(ReadPositions(&read_data))
    {
        // Skipping pre equilibrium data:
        if(containsPreEquilibriumData && read_data.timeStep < equilibrationSteps)
            continue;

        for(int i = 0; i < totalMonomers; i++)
        {
            double radialDistance = GetRadialDistance(read_data.bead_positions[i]);
            IncrementMonomerCounts(&dist_data, radialDistance);
        }
    }
    CloseReadFile(&read_data);
    FreeAllocatedMemory(&read_data);

    // Setting up the write file:
    char* writeFilePath;
    SetWriteFilePath(&writeFilePath, false, 0);
    FILE* writeFilePointer = fopen(writeFilePath, "w");
    if(writeFilePointer == NULL)
    {
        printf("The file at %s could not be opened! Terminating.\n", writeFilePath);
        exit(1);
    }

    int arraySize = dist_data.numberOfBins; // The size of the distribution array
    double distribution[arraySize];
    
    CalculateRadialDistribution(dist_data, arraySize, distribution);
    PrintDistribution(&dist_data, distribution, false, writeFilePointer); // not truncating distribution
    fclose(writeFilePointer);
    free(writeFilePath);
    FreeMonomerDistribution(&dist_data);
}

// Computes and prints the radial distribution for different regions of the polymer
void ComputeRegionalRadialDistribution(void)
{
    // Setting up reading machinery:
    char* readFilePath;
    SetReadFilePath(&readFilePath);
    InitializeSimulationData(&read_data, totalMonomers, readFilePath);
    free(readFilePath);

    // Setting up the regions struct:
    region_data regData;
    InitializeRegion(&regData, numberOfPolymers, numberOfMonomers);

    // Setting up the monomer distributions struct:
    monomer_distribution_data dist_array[regData.totalNumberOfRegions];
    for(int i = 0; i < regData.totalNumberOfRegions; i++)
    {
        int numberOfMonomersInRegion = GetNumberOfMonomersInRegion(regData, i);
        InitializeDistributionDataWithAxisLengthBinWidth(&(dist_array[i]), numberOfMonomersInRegion, radius, false, binWidth);
        // InitializeDistributionDataWithAxisLengthNumberOfBins(&dist_array[i], totalMonomers, radius, false, numberOfBins);
    }

    // Reading data and adding to distribution:
    while(ReadPositions(&read_data))
    {
        // Skipping pre equilibrium data:
        if(containsPreEquilibriumData && read_data.timeStep < equilibrationSteps)
            continue;

        for(int i = 0; i < totalMonomers; i++)
        {
            double radialDistance = GetRadialDistance(read_data.bead_positions[i]);
            int regionIndex = GetRegionID(regData, i+1); // monomer IDs are 1-indexed
            IncrementMonomerCounts(&(dist_array[regionIndex]), radialDistance);
        }
    }
    CloseReadFile(&read_data);
    FreeAllocatedMemory(&read_data);

    // Calculating and printing the distribution for each region:
    for(int i = 0; i < regData.totalNumberOfRegions; i++)
    {
        // Setting up the write files:
        char* writeFilePath;
        SetWriteFilePath(&writeFilePath, true, i);
        FILE* writeFilePointer = fopen(writeFilePath, "w");
        if(writeFilePointer == NULL)
        {
            printf("The file at %s could not be opened! Terminating.\n", writeFilePath);
            exit(1);
        }
        free(writeFilePath);
        // Calculating distribution:
        int arraySize = dist_array[i].numberOfBins; // The size of the distribution array
        double distribution[arraySize];
        
        CalculateRadialDistribution(dist_array[i], arraySize, distribution);
        PrintDistribution(&(dist_array[i]), distribution, false, writeFilePointer); // not truncating distribution
        fclose(writeFilePointer);
        FreeMonomerDistribution(&(dist_array[i]));
    }
    // Freeing the region data:
    FreeRegionMemory(&regData);

}

int main(int argc, char** argv)
{
    SetConstants(argc, argv);
    ComputeRadialDistribution();
    ComputeRegionalRadialDistribution();
    free(directory);
    printf("The radial distribution for %s Run %i has been printed.\n", architecture, runIndex);
}