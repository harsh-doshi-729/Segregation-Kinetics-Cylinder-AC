// This script calculates the monomer density of monomers along the z direction about the centre of mass of the region the monomers are a part of.

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
double axisLength; // The length of the cylinder along the z axis
bool isLengthInfinite; // A flag to indicate whether the length of the cylinder is very long as compared to the extent of the polymers
bool isBoxCentredAtOrigin = true; // A flag to indicate whether the box is centred at the origin or not
// monomer_distribution_data dist_data;

// Region information:
region_data regData; // The region data to store the region information

// File I/O information:
char* directory; // The path to the folder where the desired file is stored
simulation_data read_data;

// Directory labels:
int labelArrayLength = 2;
char* acceptedDirectoryLabels[] = {"new_segregation", "Create_Initial_States"};
bool forInitialization = false; // by default, the distribution is calculated in Create_Initial_States
int numberOfDirections = 3;
char* directions[] = {"x", "y", "z"};
int chosenDirection = 2; // z is the default direction

// Sets the directory path based on the directory label passed
void SetBaseDirectoryPath(char* label)
{
    int index = SearchStringArray(label, labelArrayLength, acceptedDirectoryLabels);
    if(index == -1)
    {
        printf("The label %s passed as the directory label could not be recognized! Please pass one of the following accepted labels: ", label);
        PrintArray(labelArrayLength, acceptedDirectoryLabels);
        exit(1);
    }
    else if(index == 0) // new_segregation passed
    {
        SetDirectoryPath(&directory, NEW_SEGREGATION, numberOfMonomers);
        forInitialization = false; // setting the flag to false for new_segregation
    }
    else if(index == 1) // Create_Initial_States passed
    {
        SetDirectoryPath(&directory, CREATE_INITIAL_STATES, numberOfMonomers);
        forInitialization = true; // setting the flag to true for Create_Initial_States
    }
}

// Reads the arguments from the command line and sets the corresponding constants:
void SetConstants(int argc, char** argv)
{
    int numberOfMandatoryArguments = 5;
    if(argc < numberOfMandatoryArguments + 1)
    {
        printf("Not enough arguments passed! Please pass the number of monomers (integer), architecture(string), the runNumber(integer), the binWidth (float), and the direction letter as arguments from the command line in the following format:\n");
        printf("./local.out <numberOfMonomers> <architecture> <runIndex> <binWidth> <direction-letter>\n");
        printf("Accepted direction letters: ");
        PrintArray(numberOfDirections, directions);
        printf("Optional argument(s) for the cylinder axis length and/or the directory label can be passed to calculate the distribution with a particular axis length in that directory.\n");
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
        char* direction = argv[5];

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
        chosenDirection = SearchStringArray(direction, numberOfDirections, directions);
        if(chosenDirection == -1)
        {
            printf("The passed argument '%s' could not be recognized as a direction. Please pass one of the following letters as the direciton argument: ", direction);
            PrintArray(numberOfDirections, directions);
            exit(1);
        }

        totalMonomers = numberOfMonomers * numberOfPolymers;
        // numberOfBins = 100;
        // binWidth = 0.1;

        containsPreEquilibriumData = false;
        // equilibrationSteps = 2 * pow(10, 6);

        // No axis length argument passed; assuming the cylinder to be infinite in length; default unless optional argument passed
        isLengthInfinite = true; // The axis length is infinite
        axisLength = GetInfiniteCylinderAxisLength(numberOfMonomers);

        // Checking for optional argument:
        if(argc > numberOfMandatoryArguments + 1)
        {
            if(argc == numberOfMandatoryArguments + 2) // Only one optional argument passed
            {
                double optionalArg = atof(argv[numberOfMandatoryArguments+1]);
                if(optionalArg == 0)
                {
                    // The optional argument passed is not a valid integer; must be the directory label
                    SetBaseDirectoryPath(argv[numberOfMandatoryArguments+1]);
                }
                else
                {
                    // A valid number passed for axisLength
                    axisLength = optionalArg;  
                    isLengthInfinite = false; // The axis length is not infinite
                }
            }
            else if(argc >= numberOfMandatoryArguments + 3) // Two or more optional arguments passed
            {
                axisLength = atof(argv[numberOfMandatoryArguments+1]);
                isLengthInfinite = false; // The axis length is not infinite
                char* label = argv[numberOfMandatoryArguments+2];
                SetBaseDirectoryPath(label); // Setting the directory based on the label passed
            }
        }
        else
            SetDirectoryPath(&directory, NEW_SEGREGATION, numberOfMonomers); // Default value of the directory path
    }
}

// Sets the read file path to the position dump file to be read from:
// Args:
// - filePathPointer: A pointer to a string where the read file path will be stored
void SetReadFilePath(char** filePathPointer)
{
    char* fileName = "visual.dump"; // default filename outside Create_Initial_States
    if(forInitialization)
        fileName = "distribution_positions.dump";
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    Append(filePathPointer, directory, folder, fileName);
}

// Sets the write file path where the radial distribution will be written:
// Args:
// - filePathPointer: A pointer to a string where the write file path will be stored
// - regionID: The ID of the region for which the distribution is being calculated (0 indexed)
void SetWriteFilePath(char** filePathPointer, int regionID)
{
    char* fileName;
    char* folderPrefix = ""; // default folder prefix for new_segregation mode
    if(forInitialization)
        folderPrefix = "monomer_distribution/"; // the preceeding folder present in Create_Initial_States run folders
    int bytes = asprintf(&fileName, "%s%s_local_monomer_distribution_reg%i.csv", folderPrefix, directions[chosenDirection], regionID + 1);
    if(bytes == -1)
    {
        printf("Memory could not be allocated for the file name!\n");
        exit(1);
    }
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    Append(filePathPointer, directory, folder, fileName);
    // Freeing the allocated memory for strings:
    free(fileName);
    free(folder);
}

// Calculates the center of mass of the region for the passed simulation data:
// Args:
// - centerOfMass: A pointer to an 2D (numberOfRegions * 3) array where the centers of mass will be stored for each region
// - read_data: The simulation data struct containing the positions of the monomers for a given timestep
// - regData: The region data struct containing the region information
// - regionID: The ID of the region for which the center of mass is being calculated
void CalculateRegionCentersOfMass(double regionCenterOfMass[][3], simulation_data read_data, region_data regData)
{
    // Initializing the center of mass to zero:
    for(int i = 0; i < regData.totalNumberOfRegions; i++)
        for(int j = 0; j < 3; j++)
            regionCenterOfMass[i][j] = 0;

    // Setting up the counters for each region:
    int counters[regData.totalNumberOfRegions]; 
    for(int i = 0; i < regData.totalNumberOfRegions; i++)
        counters[i] = 0;

    // Summing up the positions of the monomers in the region:
    for(int i = 0; i < read_data.numberOfMonomers; i++)
    {
        int regionID = GetRegionID(regData, i + 1); // regionID is 0 based
        for(int j = 0; j < 3; j++)
        {
            regionCenterOfMass[regionID][j] += read_data.bead_positions[i][j];
        }
        counters[regionID]++;
    }

    // Averaging the positions to get the center of mass:
    for(int i = 0; i < regData.totalNumberOfRegions; i++)
    {
        int numberOfMonomersInRegion = GetNumberOfMonomersInRegion(regData, i);
        // Sanity check for the number of monomers in the region:
        if(counters[i] != numberOfMonomersInRegion)
        {
            printf("The number of monomers in region %i is not equal to the expected number of monomers %i! Found: %i. Please check the region bounds in the config file.\n", 
                   i, numberOfMonomersInRegion, counters[i]);
            exit(1);
        }
        for(int j = 0; j < 3; j++)
        {
            regionCenterOfMass[i][j] /= counters[i];
        }
    }
}

// Returns the relative displacement between the two position vectors along a particular direction
// Args:
// - pos1: The first position vector (3D array)
// - pos2: The second position vector (3D array)
// - direction: The direction along which the relative distance is to be calculated (0 for x, 1 for y, 2 for z)
// Returns:
// - The relative displacement between the two position vectors along the specified direction: pos1 - pos2
double GetRelativeDisplacement(double* pos1, double* pos2, int direction)
{
    return pos1[direction] - pos2[direction]; // returns the relative distance along the specified direction
}

// Adds the relative monomer distances to the local monomer distribution for each region at a particular timestep:
// Args:
// - dist_array: An array of monomer_distribution_data structs for each region
// - read_data: The simulation data struct containing the positions of the monomers for a given timestep
// - regionCentersOfMass: A 2D (shape = numberOfregions * 3) array containing the centers of mass for each region
void AddRelativeMonomerDistancesToDistribution(monomer_distribution_data* dist_array, simulation_data read_data, double regionCentersOfMass[][3])
{
    for(int i = 0; i < read_data.numberOfMonomers; i++)
    {
        int regionID = GetRegionID(regData, i + 1); // regionID is 0 based
        double relativePosition = GetRelativeDisplacement(read_data.bead_positions[i], regionCentersOfMass[regionID], chosenDirection); // z position relative to the center of mass of the region
        IncrementMonomerCounts(&(dist_array[regionID]), relativePosition);
    }
}

// Computes and prints the monomer distribution about the center of mass of each region:
void ComputeAndPrintLocalDistribution(void)
{
    // Initializing the region struct:
    InitializeRegion(&regData, numberOfPolymers, numberOfMonomers);
    
    // Setting up the monomer distribution structs for each region:
    monomer_distribution_data* dist_array = (monomer_distribution_data*)malloc(regData.totalNumberOfRegions * sizeof(*dist_array));
    assert(dist_array);
    for(int n = 0; n < regData.totalNumberOfRegions; n++)
    {
        InitializeDistributionDataWithAxisLengthBinWidth(&(dist_array[n]), GetNumberOfMonomersInRegion(regData, n), axisLength, isBoxCentredAtOrigin, binWidth);
    }

    // Initializing the simulation data struct:
    char* readFilePath;
    SetReadFilePath(&readFilePath);
    InitializeSimulationData(&read_data, totalMonomers, readFilePath);
    free(readFilePath);

    // Initializing the centers of mass array:
    double regionCentersOfMass[regData.totalNumberOfRegions][3];

    // Reading positions and calculating distribution:
    while(ReadPositions(&read_data))
    {
        // Calculating the centers of mass for the regions:
        CalculateRegionCentersOfMass(regionCentersOfMass, read_data, regData);
        AddRelativeMonomerDistancesToDistribution(dist_array, read_data, regionCentersOfMass);
    }

    // Computing and printing the distributions:
    for(int i = 0; i < regData.totalNumberOfRegions; i++)
    {
        double distribution[dist_array[i].numberOfBins];
        CalculateDistribution(&(dist_array[i]), distribution, dist_array[i].numberOfBins);
        // Minimal printing:
        char* writeFilePath;
        SetWriteFilePath(&writeFilePath, i);
        FILE* writeFilePointer = fopen(writeFilePath, "w");
        if(writeFilePointer == NULL)
        {
            printf("The file %s could not be opened for writing! Terminating.\n", writeFilePath);
            exit(1);
        }
        PrintDistribution(&(dist_array[i]), distribution, true, writeFilePointer);
        fclose(writeFilePointer);
        free(writeFilePath); // Freeing the allocated memory for the write file path
    }

    printf("The local monomer distribution along the %s axis about the center of mass of each of the %i region(s) has been printed for %s Run %i.\n", directions[chosenDirection], regData.totalNumberOfRegions, architecture, runIndex);

    // Freeing the allocated memory for other structs:
    FreeAllocatedMemory(&read_data);
    FreeRegionMemory(&regData);
    for(int i = 0; i < regData.totalNumberOfRegions; i++)
        FreeMonomerDistribution(&(dist_array[i]));
    free(dist_array); // Freeing the allocated memory for the distribution array
}

// Debugging function: prints the centers of mass of the regions to the standard output
void PrintRegionCentersOfMass(double regionCentersOfMass[][3], int totalRegions)
{
    printf("Region Centers of Mass:\n");
    for(int i = 0; i < totalRegions; i++)
    {
        printf("Region %i: (%.2lf, %.2lf, %.2lf)\n", i + 1, regionCentersOfMass[i][0], regionCentersOfMass[i][1], regionCentersOfMass[i][2]);
    }
    printf("\n");
}

// Debugging function: printing the calculated relative distances to verify the calculation
void PrintRelativeDistances(FILE* writeFilePointer, simulation_data read_data, double regionCentersOfMass[][3])
{
    for(int i = 0; i < read_data.numberOfMonomers; i++)
    {
        int regionID = GetRegionID(regData, i + 1); // regionID is 0 based
        double relativePosition = GetRelativeDisplacement(read_data.bead_positions[i], regionCentersOfMass[regionID], chosenDirection); // position relative to the center of mass of the region
        fprintf(writeFilePointer, "%i,%lf\n", i + 1, relativePosition);
    }
}

int test_main(int argc, char** argv)
{
    // Setting the constants from the command line arguments:
    SetConstants(argc, argv);

    // Initializing the region data:
    InitializeRegion(&regData, numberOfPolymers, numberOfMonomers);

    // Initializing the simulation data:
    char* readFilePath;
    SetReadFilePath(&readFilePath);
    InitializeSimulationData(&read_data, totalMonomers, readFilePath);
    free(readFilePath);

    // Printing the constants for debugging:
    printf("Number of monomers: %i\n", numberOfMonomers);
    printf("Architecture: %s\n", architecture);
    printf("Run index: %i\n", runIndex);
    printf("Bin width: %lf\n", binWidth);
    printf("Axis length: %.2lf\n", axisLength);
    printf("Is length infinite: %s\n", isLengthInfinite ? "true" : "false");
    printf("Directory: %s\n", directory);
    printf("Chosen Direction: %s\n", directions[chosenDirection]);

    // Computing and printing the local distribution:
    // ComputeAndPrintLocalDistribution();
    // Debugging: printing the relative distances
    char* writeFilePath;
    SetWriteFilePath(&writeFilePath, 0); // passing 0 as the region ID for the write file path
    FILE* writeFilePointer = fopen(writeFilePath, "w");
    if(writeFilePointer == NULL)
    {
        printf("The file %s could not be opened for writing! Terminating.\n", writeFilePath);
        exit(1);
    }
    double regionCentersOfMass[regData.totalNumberOfRegions][3];
    
    ReadPositions(&read_data); // Reading the positions for the first timestep
    CalculateRegionCentersOfMass(regionCentersOfMass, read_data, regData);
    PrintRegionCentersOfMass(regionCentersOfMass, regData.totalNumberOfRegions); // Debugging: printing the centers of mass
    PrintRelativeDistances(writeFilePointer, read_data, regionCentersOfMass);

    return 0;
}

int main(int argc, char** argv)
{
    SetConstants(argc, argv);
    ComputeAndPrintLocalDistribution();
}
