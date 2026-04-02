// This script is used to calculate the mean squared displacement (MSD) of the COMs of the polymers while segregating
// The hope is that the MSD will change slope after the polymers have achieved segregation, and thus, can give an indication of the segregation time

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

// Global variables:
int numberOfPolymers = 2;
int numberOfMonomers; // The number of monomomers present in each polymer

char* architecture; // The name of the architecture being used
int runIndex; // The index of the run whose data is to be processed

char* directory; // The directory where the runs are stored

// Timing information:
int timeInterval; // The time interval between successive position snapshots in the COM file
int totalTimeSteps; // The total number of timesteps in the simulation; assuming the run starts from 0
int arrayLength; // The length of the array to store the COM values at each timestep

// Mode information:
int numberOfAcceptedModes = 3;
char* acceptedModes[] = {"COM", "relative_COM", "z_COM"};
int chosenMode; // The index of the chosen mode

// Sets the variables for the timing information based on the number of monomers
void SetTimingConstants(int numberOfMonomers)
{
    switch(numberOfMonomers)
    {
        case 200:
            timeInterval = 1000;
            totalTimeSteps = 4 * pow(10, 7);
            break;
        case 500:
            timeInterval = 1000;
            totalTimeSteps = 2.5 * pow(10, 8);
            break;
        default:
            printf("ERROR: The number of monomers %i is not recognized! Accepted values are 200 and 500.\n", numberOfMonomers);
            exit(1);
    }
}


// Setting the constants from command line arguments
void SetConstants(int argc, char** argv)
{
    int numberOfMandatoryArguments = 4; 
    if(argc < numberOfMandatoryArguments + 1) // +1 to include the program name
    {
        printf("Not enough arguments passed! Please pass the number of monomers (integer), the architecture name (string), the run index (integer), and the operation mode (string) as arguments to the command line.\n");
        printf("Accepted modes: ");
        PrintArray(numberOfAcceptedModes, acceptedModes);
        printf("For example: ./a.out 200 Arc2 1 COM\n");
        exit(1);
    }
    else
    {
        // Reading arguments:
        numberOfMonomers = atoi(argv[1]);
        if(numberOfMonomers == 0)
        {
            printf("ERROR: The argument %s could not be converted to a valid integer for the number of monomers! Please pass a positive integer.\n", argv[1]);
            exit(1);
        }
        architecture = argv[2];
        runIndex = atoi(argv[3]);
        if(runIndex <= 0)
        {
            printf("ERROR: The argument %s could not be converted to a valid positive integer for the run index! Please pass a positive integer.\n", argv[3]);
            exit(1);
        }

        char* mode = argv[4];
        chosenMode = SearchStringArray(mode, numberOfAcceptedModes, acceptedModes);
        if(chosenMode == -1)
        {
            printf("ERROR: The passed mode '%s' was not recognized. Please pass one of the following: ", mode);
            PrintArray(numberOfAcceptedModes, acceptedModes);
            exit(1);
        }

        // Setting other constants:
        SetDirectoryPath(&directory, NEW_SEGREGATION, numberOfMonomers);
        // Setting timing constants:
        SetTimingConstants(numberOfMonomers);
        arrayLength = totalTimeSteps / timeInterval; 
    }
}

// Sets the path to the file where the MSD data will be written
void SetMSDWriteFilePath(char** filePathPointer, char* directory, char* architecture, int runIndex)
{
    char* fileName;
    int bytes = asprintf(&fileName, "%s_MSD.csv", acceptedModes[chosenMode]);
    if(bytes == -1)
    {
        printf("ERROR: Memory could not be allocated for the MSD file name.\n");
        exit(1);
    }
    SetCustomWriteFilePath(filePathPointer, directory, architecture, runIndex, fileName);
    free(fileName);
}

// Returns the squared displacement between two position vectors
double CalculateSquaredDisplacement(double position1[3], double position2[3])
{
    double squaredDisplacement = 0.0;
    for(int i = 0; i < 3; i++)
    {
        squaredDisplacement += pow(position1[i] - position2[i], 2);
    }
    return squaredDisplacement;
}

// Calculates the squared COM displacement along a single axis/direction between two position vectors
// Args:
// - position1: The first position vector
// - position2: The second position vector
// - axis: The axis along which the squared displacement is to be calculated; 0 for x, 1 for y, 2 for z
double CalculateSquaredDisplacementAlongAxis(double position1[3], double position2[3], int axis)
{
    return pow(position1[axis] - position2[axis], 2);
}

// Calculate relative COM displacement between the two polymers
// Args:
// - relative_COM_Array: A 2D array to store the relative COM displacement vectors for each timestep
// - COM_Array: A 3D double array that contains the COM position vectors at each timestep for each polymer
void CalculateRelativeCOM(double relative_COM_Array[arrayLength][3], double COM_Array[numberOfPolymers][arrayLength][3])
{
    for(int i = 0; i < arrayLength; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            relative_COM_Array[i][j] = COM_Array[0][i][j] - COM_Array[1][i][j];
        }
    }
}

// Reads the COM data from the specified file and stores it in the provided array
// Args:
// - filePath: The path to the file containing the COM data
// - arrayLength: The length of the array to store the COM data; the number of timesteps at which COM data is recorded
// - COM_data: A 2D array to store the COM data; should be of size arrayLength x 3
void ReadCOMData(char* filePath, int arrayLength, double COM_data[][3])
{
    FILE* filePointer = fopen(filePath, "r");
    if(filePointer == NULL)
    {
        printf("ERROR: The file %s could not be opened for reading!\n", filePath);
        exit(1);
    }

    // Reading the data:
    int maxLineSize = 100;
    char buffer[maxLineSize];
    // Skipping header line:
    fgets(buffer, maxLineSize, filePointer);
    // Counting number of lines read:
    int counter = 0;
    while(true)
    {
        fgets(buffer, maxLineSize, filePointer);

        if(feof(filePointer))
        {
            printf("COM file EOF reached.\n");
            break;
        }

        int successes = sscanf(buffer, "%*i %lf %lf %lf\n", &COM_data[counter][0], &COM_data[counter][1], &COM_data[counter][2]);
        if(successes != 3)
        {
            printf("ERROR: The line %s could not be parsed correctly! Only %i values were read instead of 3.\n", buffer, successes);
            exit(1);
        }
        counter++;
        // Sanity check:
        if(counter > arrayLength)
        {
            printf("ERROR: The number of lines read in the COM file has exceeded the expected number of lines %i!\n", arrayLength);
            exit(1);
        }
    }
    fclose(filePointer);
}

// Calculates the mean squared displacement (MSD) of the polymers based on their COM data, and writes it to a file
// Args:
// - arrayLength: The length of the COM_data
// - COM_data: A 2D array containing the COM data of the polymers; should be of size arrayLength x 3
// - writeFilePath: The path to the file where the MSD data should be written
// - timeInterval: The time interval between successive COM data points
void CalculateAndWriteMSD(void)
{
    // Setting up write file:
    char* writeFilePath;
    SetMSDWriteFilePath(&writeFilePath, directory, architecture, runIndex);
    // Opening file:
    FILE* filePointer = fopen(writeFilePath, "w");
    if(filePointer == NULL)
    {
        printf("ERROR: The file %s could not be opened for writing!\n", writeFilePath);
        exit(1);
    }
    // Writing header:
    fprintf(filePointer, "Delta Time, MSD (%s)\n", acceptedModes[chosenMode]);

    // Reading the COM timeseries of the two polymers:
    double COM_Array[numberOfPolymers][arrayLength][3]; // 3D array to store the COM data of both polymers
    for(int i = 0; i < numberOfPolymers; i++)
    {
        char* COM_filePath;
        char* fileName;
        int bytes = asprintf(&fileName, "com%i.dat", i+1);
        if(bytes == -1)
        {
            printf("ERROR: Memory could not be allocated for the COM file name.\n");
            exit(1);
        }
        SetCustomWriteFilePath(&COM_filePath, directory, architecture, runIndex, fileName);
        free(fileName);
        ReadCOMData(COM_filePath, arrayLength, COM_Array[i]);
        free(COM_filePath);
    }
    // Debugging: Printing squared displacement for all windows of a particular length
    int numberOfWindowLengths = 6;
    int windowLengths[] = {50000 / timeInterval, 100000 / timeInterval, 800000 / timeInterval, 1000000 / timeInterval, 2000000 / timeInterval, 8000000 / timeInterval}; // in time index
    FILE* debugFilePointers[numberOfWindowLengths];
    char* debugFileNames[] = {"relative_COM_MSD_0005_window.csv", "relative_COM_MSD_001_window.csv", "relative_COM_MSD_008_window.csv", "relative_COM_MSD_01_window.csv", "relative_COM_MSD_02_window.csv", "relative_COM_MSD_08_window.csv"};
    bool isWindowLength[numberOfWindowLengths];
    for(int i = 0; i < numberOfWindowLengths; i++)
    {
        char* debugFilePath;
        SetCustomWriteFilePath(&debugFilePath, directory, architecture, runIndex, debugFileNames[i]);
        debugFilePointers[i] = fopen(debugFilePath, "w");
        free(debugFilePath);
        fprintf(debugFilePointers[i], "Start Time, Squared Displacement\n");
    }
    // Calculating the relative COM array if needed:
    double relativeCOM[arrayLength][3];
    if(acceptedModes[chosenMode] == "relative_COM")
        CalculateRelativeCOM(relativeCOM, COM_Array);
    // Calculating the MSD for different delta times:
    for(int deltaIndex = 0; deltaIndex < arrayLength; deltaIndex++)
    {
        int deltaTime = deltaIndex * timeInterval;
        double MSD = 0.0;
        int count = 0; // The number of valid terms added to the MSD sum
        // Debugging: Setting window lengths:
        for(int i = 0; i < numberOfWindowLengths; i++)
        {
            if(deltaIndex == windowLengths[i])
                isWindowLength[i] = true;
            else
                isWindowLength[i] = false;
        }
        for(int startIndex = 0; startIndex + deltaIndex < arrayLength; startIndex++)
        {
            if(acceptedModes[chosenMode] == "COM")
            {
                // Calculating the squared displacement for both polymers:
                for(int j = 0; j < numberOfPolymers; j++)
                {
                    double squaredDisplacement = CalculateSquaredDisplacement(COM_Array[j][startIndex], COM_Array[j][startIndex + deltaIndex]);
                    MSD += squaredDisplacement;
                    count++;
                }
            }
            else if(acceptedModes[chosenMode] == "relative_COM")
            {
                // Calculating the squared displacement of the COMs of the two polymers
                double squaredDisplacement = CalculateSquaredDisplacement(relativeCOM[startIndex], relativeCOM[startIndex + deltaIndex]);
                MSD += squaredDisplacement;
                for(int k = 0; k < numberOfWindowLengths; k++)
                {
                    if(isWindowLength[k])
                        fprintf(debugFilePointers[k], "%i, %.6lf\n", startIndex * timeInterval, squaredDisplacement);
                }
                count++;
            }
            else if(acceptedModes[chosenMode] == "z_COM")
            {
                // Calculating the squared displacement along the z axis for both polymers:
                for(int j = 0; j < numberOfPolymers; j++)
                {
                    double squaredDisplacement = CalculateSquaredDisplacementAlongAxis(COM_Array[j][startIndex], COM_Array[j][startIndex + deltaIndex], 2); // 2 for z axis
                    MSD += squaredDisplacement;
                    count++;
                }
            }
            else
            {
                printf("ERROR: The chosen mode '%s' was not recognized!\n", acceptedModes[chosenMode]);
                exit(1);
            }
        }
        if(count > 0)
            MSD /= count; // Taking the average
        else // If no valid terms were found, set MSD to 0
            MSD = 0.0;

        // Writing to file:
        fprintf(filePointer, "%i, %.6lf\n", deltaTime, MSD);
    }
    fclose(filePointer);
    // Closing debugging files:
    for(int i = 0; i < numberOfWindowLengths; i++)
        fclose(debugFilePointers[i]);
    free(writeFilePath);
    // Printing success message:
    printf("%s MSD calculation completed successfully for %i polymer(s) of %s Run %i.\n", acceptedModes[chosenMode], numberOfPolymers, architecture, runIndex);
}

int main(int argc, char** argv)
{
    SetConstants(argc, argv);
    
    // timing the process:
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);
    printf("Start (Number of Monomers = %i, Architecture = %s, Run Index = %i, Special Simulation = %s): %s\n", numberOfMonomers, architecture, runIndex, SPECIAL_SIMULATION, asctime(tm));
    clock_t start = clock();

    CalculateAndWriteMSD();

    //printing timing information:
    t = time(NULL);
    tm = localtime(&t);
    clock_t finish = clock();
    long int elapsedTime = (double)(finish - start)/ CLOCKS_PER_SEC;
    printf("Finish: %s\nElapsed (CPU) Time: %lih %lim %lis\n", asctime(tm), elapsedTime/3600, (elapsedTime/60) % 60, elapsedTime%60);

    return 0;
}