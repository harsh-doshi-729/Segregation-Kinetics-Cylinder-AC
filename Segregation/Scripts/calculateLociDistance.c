// This script is used to calculate the distance between two particular monomers (loci) present in a simulation
// The distance is calculated at each timestep that is dumped

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

// Global variables:
int numberOfMonomers; // The number of monomers in a single polymer
int numberOfPolymers = 2; // The number of polymers in the system
int totalMonomers; // The total number of monomers in the system
char* architecture;
int runIndex;

// File information:
int numberOfKeywords = 3;
char* acceptedDirectoryKeywords[3] = {"new_segregation", "create_initial_states", "custom"};
char* keywordInstructions[3] = {"1 additional argument; the run index (integer) for which the computation should be performed",
                                "1 additional argument; the run index (integer) for which the computation should be performed",
                                "1 additional argument; the custom folder path (string) where the desired dump file exists"};
char* chosenKeyword;
bool customDirectory; // A flag to indicate whether the custom keyword is chosen
char* directory; // The path to the directory where all the concatenation dump files are stored
char* fileName;

// Distance information:
int indexPair[2];

// Sets the directory variable based on the directory keyword passed
void SetDirectory(int numberOfMandatoryArguments, char* chosenKeyword, char** argv)
{
    if(strcasecmp(chosenKeyword, acceptedDirectoryKeywords[numberOfKeywords - 1]) == 0) // custom
    {
        char* customPath = argv[numberOfMandatoryArguments+1];
        Append(&directory, customPath, "", "");
        customDirectory = true;
    }
    else // new_segregation and create_initial_states
    {
        customDirectory = false;
        int runIndex = atoi(argv[numberOfMandatoryArguments+1]);
        char* folderPath;
        char* pathPrefix;
        if(strcasecmp(chosenKeyword, "new_segregation") == 0)
            pathPrefix = NEW_SEGREGATION;
        else if(strcasecmp(chosenKeyword, "create_initial_states") == 0)
            pathPrefix = CREATE_INITIAL_STATES;
        SetDirectoryPath(&folderPath, pathPrefix, numberOfMonomers);
        char* runLabel;
        int bytes = asprintf(&runLabel, "/run%i/", runIndex);
        if(bytes == -1)
        {
            printf("ERROR: Memory could not be allocated for the run label!\n");
            exit(1);
        }
        Append(&directory, folderPath, architecture, runLabel);
        free(folderPath);
        free(runLabel);
    }
}

// Sets the relevant constants for this script after reading the inputs from the command line
void SetConstants(int argc, char** argv)
{
    int numberOfMandatoryArguments = 6; // atleast these arguments excluding keyword arguments
    char* directory;
    if(argc < numberOfMandatoryArguments + 1) // not enough arguments were passed
    {
        printf("Not enough arguments passed! Please pass the numberOfMonomers(integer), the architecture(string), the pair of monomer indices between which distance should be calculated, the dump file name (string), the directory keyword (string), and any subsequent required keywords as arguments from the command line.\n");
        printf("The directory keyword must be chosen from the following: ");
        PrintArray(numberOfKeywords, acceptedDirectoryKeywords);
        PrintKeywordInstructions(numberOfKeywords, acceptedDirectoryKeywords, keywordInstructions);
        printf("For example, ./a.out 200 Arc1_10 151 50 visual.dump new_segregation 1\n");
        exit(1);
    }
    else
    {
        numberOfMonomers = atoi(argv[1]);
        totalMonomers = numberOfPolymers * numberOfMonomers;
        architecture = argv[2];
        indexPair[0] = atoi(argv[3]);
        indexPair[1] = atoi(argv[4]);
        fileName = argv[5];
        chosenKeyword = argv[6];
        if(SearchStringArray(chosenKeyword, numberOfKeywords, acceptedDirectoryKeywords) == -1)
        {
            printf("The passed directory keyword '%s' was not recognized. Please pass a keyword from the following list: ", chosenKeyword);
            PrintArray(numberOfKeywords, acceptedDirectoryKeywords);
            PrintKeywordInstructions(numberOfKeywords, acceptedDirectoryKeywords, keywordInstructions);
            exit(1);
        }
        SetDirectory(numberOfMandatoryArguments, chosenKeyword, argv);
    }
}

// Calculate the absolute differences between the coordinates of the two monomers along each Cartesian direction
// Args:
//      differences: A double array of length 3 to store teh difference along each direction
//      read_data: The struct to read the position dump data
//      indexPair: An integer array of size 2 containing a pair of monomer indices (starting from 1) between which distance is to be calculated
void CalculateCoordinateDifference(double* differences, simulation_data read_data, int indexPair[2])
{
    for(int i = 0; i < 3; i++) // iterating over spatial dimensions
    {
        differences[i] = fabs(read_data.bead_positions[indexPair[0]-1][i] - read_data.bead_positions[indexPair[1]-1][i]);
    }
}

// Calculates the magnitude of the three dimensional vector passed
double CalculateMagnitude(double* vector)
{
    double distance = 0;
    for(int i = 0; i < 3; i++) // iterating over spatial dimensions
        distance += pow(vector[i], 2);
    distance = sqrt(distance);
    return distance;
}

// Reads the position dump file and prints the distances between the two indexed monomers for all timesteps
// Args:
//      totalNumberOfMonomers: the total number of monomers expected to be present in the system
//      indexPair: an integer array of size 2 containing a pair of monomer indices (starting from 1) between which distance is to be calculated
//      readFilePath: file path to the position dump file of the simulation
//      writeFilePath: file path to the file to which the distances are to be written
void ReadAndCalculateDistance(int totalNumberOfMonomers, int indexPair[2], char* readFilePath, char* writeFilePath)
{
    // Setting up read file:
    simulation_data read_data; // The struct for the LAMMPS Position dump file
    InitializeSimulationData(&read_data, totalNumberOfMonomers, readFilePath);

    // Setting up write file:
    FILE* writeFilePointer = fopen(writeFilePath, "w");
    if(writeFilePointer == NULL)
    {
        printf("ERROR: The write file %s could not be opened!\n", writeFilePath);
        exit(1);
    }
    fprintf(writeFilePointer, "Timestep, x-diff, y-diff, z-diff, distance\n"); // Writing header

    // Reading the file and calculating distance, their statistics:
    double differences[3];
    double means[4] = {0, 0, 0, 0}; 
    double stds[4] = {0, 0, 0, 0};
    int counter = 0;

    while(ReadPositions(&read_data))
    {
        CalculateCoordinateDifference(differences, read_data, indexPair);
        // Writing data to file:
        fprintf(writeFilePointer, "%i, ", read_data.timeStep);
        for(int i = 0; i < 3; i++) // iterating over spatial dimensions
            fprintf(writeFilePointer, "%.6lf, ", differences[i]);
        double distance = CalculateMagnitude(differences);
        fprintf(writeFilePointer, "%.6lf\n", distance);
        // Calculating statistics:
        for(int i = 0 ; i < 3; i++)
        {
            means[i] += differences[i];
            stds[i] += pow(differences[i], 2);
        }
        means[3] += distance;
        stds[3] += pow(distance, 2);
        counter++;
    }

    // Calculating and Printing statistics:
    printf("\nPrinting statistics for distance between monomers %i and %i for %s(%i):\n", indexPair[0], indexPair[1], architecture, numberOfMonomers);
    printf("Headers: x-diff, y-diff, z-diff, Distance\n");
    for(int i = 0; i < 4; i++)
    {
        means[i] = means[i] / counter;
        stds[i] = sqrt(stds[i] / (counter-1) - counter / (counter - 1) * pow(means[i], 2)); // using the N-1 degrees of freedom
    }
    // Printing means:
    bothprintf(writeFilePointer, "Mean: ");
    for(int i = 0; i < 3; i++)
        bothprintf(writeFilePointer, "%.6lf, ", means[i]);
    bothprintf(writeFilePointer, "%.6lf\n", means[3]);
    // Printing STDs:
    bothprintf(writeFilePointer, "Std. Dev.: ");
    for(int i = 0; i < 3; i++)
        bothprintf(writeFilePointer, "%.6lf, ", stds[i]);
    bothprintf(writeFilePointer, "%.6lf\n", stds[3]);
    fclose(writeFilePointer);
    
    // Freeing memory:
    FreeAllocatedMemory(&read_data);
    printf("Calculated distance data between monomers %i and %i for %s(%i).\n",  indexPair[0], indexPair[1], architecture, numberOfMonomers);
}

// Sets the write file name for a particular pair of monomer indices (loci)
void SetWriteFileName(char** fileNamePointer, int indexPair[2])
{
    int bytes = asprintf(fileNamePointer, "distance-%i-%i.csv", indexPair[0], indexPair[1]);
    if(bytes == -1)
    {
        printf("ERROR: Memory could not be allocated for the distance write file name.\n");
        exit(1);
    }
}

// Executes the calculation of the distances between the monomer defined in this script
void ExecuteCalculation(void)
{
    char* readFilePath;
    Append(&readFilePath, directory, fileName, "");
    char* writeFileName;
    SetWriteFileName(&writeFileName, indexPair);
    char* writeFilePath;
    Append(&writeFilePath, directory, writeFileName, "");
    free(writeFileName);

    int numberOfPolymers = ReadNumberOfPolymers(numberOfMonomers, readFilePath);
    totalMonomers = numberOfPolymers * numberOfMonomers;
    ReadAndCalculateDistance(totalMonomers, indexPair, readFilePath, writeFilePath);
    
    free(readFilePath);
    free(writeFilePath);
}

int main(int argc, char** argv)
{
    SetConstants(argc, argv);

    // // Printing read information:
    // printf("Number of Monomers = %i\nArchitecture = %s\n", numberOfMonomers, architecture);
    // printf("Monomer Indices = (%i, %i)\n", indexPair[0], indexPair[1]);
    // printf("Chosen Keyword = %s\nDirectory = %s\nFile Name = %s\n", chosenKeyword, directory, fileName);
    // printf("Is Directory custom = %s\n", customDirectory ? "true": "false");
    // exit(0);
    
    // Timing the calculation:
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);
    printf("Start: %s\n", asctime(tm));
    clock_t start = clock();

    ExecuteCalculation();
    free(directory);

    t = time(NULL);
    tm = localtime(&t);
    clock_t finish = clock();
    long int elapsedTime = (double)(finish - start)/ CLOCKS_PER_SEC;
    printf("Finish: %s\nElapsed (CPU) Time: %lih %lim %lis\n", asctime(tm), elapsedTime/3600, (elapsedTime/60) % 60, elapsedTime%60);

    return 0;
}