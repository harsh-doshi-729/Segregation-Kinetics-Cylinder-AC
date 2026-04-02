// This script calculates the velocity (time derivative of positions) along z direction of each monomer from the LAMMPS position dump files

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
#include REGION


// Steps:
// 1. Read the data: The position dump files must be read, and the positions stored
// 2. Calculate velocity as the time derivative of position: v(t) = (x(t+dt) - x(t)) / dt; or some more complicated scheme: centered finite difference v(t) = {x(t+dt)-x(t-dt)}/2dt
// 3. Print the COM velocities as a function of time

// Global Variables:
int numberOfPolymers = 2;
int numberOfMonomers; // The number of monomomers present in each polymer
int totalMonomers; // The total number of monomers in the system
char* architecture; // The name of the architecture being used
int runIndex; // The index of the run whose data is to be processed
int numberOfRuns = 50; // The total number of runs present
int timeInterval; // The time interval between successive position snapshots in the dump file; in units of iterations

// Which velocity to calculate?
char* acceptedModes[] = {"relative-com"};
int numberOfAcceptedModes = 1;
int chosenMode; // The index of the mode to be calculated
int regionPair[2] = {1, 2}; // A pair of regions between which the displacement is calculated; default is the polymerical region pair: 1, 2
bool useConfigRegions = false; // A flag to indicate whether to use the regions in the config file; false by default

// Directory information:
char* directory;

void SetConstants(int argc, char** argv)
{
    int numberOfMandatoryArguments = 4; 
    if(argc < numberOfMandatoryArguments + 1) // +1 to include the program name
    {
        printf("Not enough arguments passed! Please pass the number of monomers (integer), the architecture name (string), the run index (integer), and the velocity mode (string) as arguments to the command line.\n");
        printf("Accepted velocity modes: ");
        PrintArray(numberOfAcceptedModes, acceptedModes);
        printf("For example: ./a.out 200 Arc2 1 relative-com\n");
        printf("Optionally, you can specify a pair of integer regions as listed in the config file (e.g., 1 2) to calculate the displacement between them.\n");
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
        totalMonomers = numberOfPolymers * numberOfMonomers;

        architecture = argv[2];
        runIndex = atoi(argv[3]);
        if(runIndex == 0)
        {
            printf("ERROR: The argument %s could not be converted to a valid integer for the run index! Please pass a positive integer.\n", argv[3]);
            exit(1);
        }

        char* modeArg = argv[4];
        chosenMode = SearchStringArray(modeArg, numberOfAcceptedModes, acceptedModes);
        if(chosenMode == -1)
        {
            printf("ERROR: The velocity mode identifier '%s' was not recognized. Please pass one of the following accepted identifiers: ", modeArg);
            PrintArray(numberOfAcceptedModes, acceptedModes);
            exit(1);
        }

        // Setting directory:
        SetDirectoryPath(&directory, NEW_SEGREGATION, numberOfMonomers);

        // Optional arguments:
        if(argc > numberOfMandatoryArguments + 1) // Atleast one optional argument passed
        {
            if(argc < numberOfMandatoryArguments + 1 + 2)
            {
                printf("ERROR: Not enough optional arguments passed! Please specify a pair of integer regions (e.g., 1 2) to calculate the displacement between them.\n");
                exit(1);
            }
            else
            {
                regionPair[0] = atoi(argv[5]);
                regionPair[1] = atoi(argv[6]);
                if(regionPair[0] == 0 || regionPair[1] == 0)
                {
                    printf("ERROR: The arguments %s and %s could not be converted to valid region indices for the region pair! Please pass positive integers.\n", argv[5], argv[6]);
                    exit(1);
                }
                useConfigRegions = true;
            }
        }
    }
}

// Sets the file path to the write file where the calculated velocities will be written
void SetWriteFilePath(char** filePathPointer, char* velocityMode, int regionPair[2])
{
    char* fileName;
    int bytes;
    if(useConfigRegions)
        bytes = asprintf(&fileName, "velocity_%s_reg_%i_%i_run%i.csv", velocityMode, regionPair[0], regionPair[1], runIndex);
    else
        bytes = asprintf(&fileName, "velocity_%s_run%i.csv", velocityMode, runIndex);
    if(bytes == -1)
    {
        printf("ERROR: Memory could not be allocated for the write file name!\n");
        exit(1);
    }
    SetCustomWriteFilePath(filePathPointer, directory, architecture, runIndex, fileName);
    free(fileName);
}

// Calculates the displacement vector as per the chosen mode and stores it in the passed vector
// Args:
// - vector: A double array of size 3 to store the displacement vector
// - chosenMode: An integer indicating the mode of displacement to be calculated
// - read_data: The simulation_data struct that contains the read monomer positions for the current timestep
// - regData: The region_data struct that stores information about the regions the system is divided in
// - region1: The index of the first region; indexed from 1
// - region2: The index of the second region; indexed from 1
void GetGenericDisplacement(double* vector, int chosenMode, simulation_data read_data, region_data regData, int region1, int region2)
{
    switch(chosenMode)
    {
        case 0: // relative-com
        {
            int numberOfRegions = regData.totalNumberOfRegions;
            double coms[numberOfRegions][3];
            for(int i = 0; i < numberOfRegions; i++)
            {
                CalculateRegionCoM(coms[i], read_data, regData, i);
            }
            for(int j = 0; j < 3; j++)
            {
                vector[j] = coms[region1-1][j] - coms[region2-1][j];
            }
            break;
        }
        default:
        {
            printf("ERROR: The chosen mode %i was not recognized! Terminating.\n", chosenMode);
            exit(1);
        }
    }
}

// Calculates the time derivative of the distance (chosen by mode) and prints the velocities as a function of time to a file
void CalculateAndPrintVelocities(int regionPair[2])
{
    // Initializing region:
    region_data regData;
    if(useConfigRegions)
        InitializeRegion(&regData, numberOfPolymers, numberOfMonomers); // Regions according to the config file
    else
        InitializeRegionManually(&regData, numberOfPolymers, numberOfMonomers, true, false); // Polymeric regions only

    // Setting file paths:
    char* readFilePath;
    SetPositionDumpFilePath(&readFilePath, directory, architecture, runIndex, false); // Only for segregation runs, for now
    char* writeFilePath;
    SetWriteFilePath(&writeFilePath, acceptedModes[chosenMode], regionPair);

    // Reading position file:
    simulation_data read_data;
    InitializeSimulationData(&read_data, totalMonomers, readFilePath);
    free(readFilePath);
    
    // Opening write file:
    FILE* writeFile = fopen(writeFilePath, "w");
    if(writeFile == NULL)
    {
        printf("ERROR: The write file %s could not be opened!\n", writeFilePath);
        exit(1);
    }
    free(writeFilePath);

    // Printing header to write file:
    fprintf(writeFile, "Time (iterations), Velocity_x, Velocity_y, Velocity_z\n");

    // Setting up a (moving) time window in which the displacement vectors can be stored:
    double firstDisplacement[3]; // distance at t-1
    double secondDisplacement[3]; // distance at t
    double thirdDisplacement[3]; // distance at t+1
    // Iterating over all timesteps:
    bool firstTimestep = true; // flag to indicate whether the first time step has just been read
    bool secondTimestep = false; // flag to indicate whether the second time step has just been read
    while(ReadPositions(&read_data))
    {
        double velocity[3];
        if(firstTimestep) // if this is the first timestep being read
        {
            GetGenericDisplacement(firstDisplacement, chosenMode, read_data, regData, regionPair[0], regionPair[1]); // Getting displacement at t-1
            firstTimestep = false;
            secondTimestep = true; // Set this flag to true for the next iteration of the while loop
            continue; // skipping the rest of the loop for the first timestep
        }
        else if(secondTimestep)
        {
            timeInterval = read_data.printingInterval; // setting the time interval between successive snapshots
            GetGenericDisplacement(secondDisplacement, chosenMode, read_data, regData, regionPair[0], regionPair[1]); // Getting displacement at t
            secondTimestep = false; // Reset this flag for the next iteration
            for (int j = 0; j < 3; j++)
            {
                // Calculate first derivative - forward derivative:
                velocity[j] = (secondDisplacement[j] - firstDisplacement[j]) / (timeInterval); // in units of distance per iteration
            }
            
        }
        else
        {
            // Because the derivative calculation is done via centred finite difference (data from three points in time is required),
            // the third displacement is the one read from the current timestep
            GetGenericDisplacement(thirdDisplacement, chosenMode, read_data, regData, regionPair[0], regionPair[1]); // Getting displacement at t+1

            // The velocity is calculated for the previously read positions i.e read_data.timestep-1; which is considered as time t
            for(int j = 0; j < 3; j++)
            {
                // Calculate derivative - centred derivative:
                velocity[j] = (thirdDisplacement[j] - firstDisplacement[j]) / (2 * timeInterval); // in units of distance per iteration
                // Shifting time window / Updating the displacements for the next timestep
                firstDisplacement[j] = secondDisplacement[j]; // updating first displacement
                secondDisplacement[j] = thirdDisplacement[j]; // updating second displacement
            }
        }
        
        // Printing velocity to file for the previous step:
        fprintf(writeFile, "%i, %lf, %lf, %lf\n", read_data.timeStep - timeInterval, velocity[0], velocity[1], velocity[2]);
    }
    // Calculating the velocity for the last step: backward difference
    double lastVelocity[3];
    for(int j = 0; j < 3; j++)
    {
        lastVelocity[j] = (secondDisplacement[j] - firstDisplacement[j]) / (timeInterval);
    }
    // Printing the last velocity:
    fprintf(writeFile, "%i, %lf, %lf, %lf\n", read_data.timeStep, lastVelocity[0], lastVelocity[1], lastVelocity[2]);
    fclose(writeFile);

    // Printing success message:
    char* regionLabel = "";
    if(useConfigRegions)
    {
        int bytes = asprintf(&regionLabel, "for regions %i and %i ", regionPair[0], regionPair[1]);
        if(bytes == -1)
        {
            printf("ERROR: Could not allocate memory for the region label.\n");
            exit(1);
        }
    }
    printf("Calculated and printed %s velocities %sfor %s Run %i.\n", acceptedModes[chosenMode], regionLabel, architecture, runIndex);
    if(useConfigRegions)
        free(regionLabel);

    // Freeing memory of structs:
    FreeAllocatedMemory(&read_data);
    FreeRegionMemory(&regData);
}

// Reads the velocity dump of a simulations and calculates the relative COM velcoties at each timestep
void ReadAndPrintVelocities(void)
{
    char* velocityFileDump = "/home/harsh/ComputationalPhysics/Polymer_Physics/LAMMPS_runs/cluster-data/new_segregation/b200/Previous_Attempts/fene_recenter/Arc1_10/run1/velocity.dump"; // setting it manually for now

    // Using the Position file struct to read the velocity:
    simulation_data read_data;
    InitializeSimulationData(&read_data, totalMonomers, velocityFileDump);
    char* writeFilePath = "/home/harsh/ComputationalPhysics/Polymer_Physics/LAMMPS_runs/cluster-data/new_segregation/b200/Previous_Attempts/fene_recenter/Arc1_10/run1/com_velocity.csv";

    FILE* writeFile = fopen(writeFilePath, "w");
    if(writeFile == NULL)
    {
        printf("ERROR: The manual write file path could not be opened.\n");
        exit(1);
    }
    // Printing header:
    fprintf(writeFile, "Timestep, v_x, v_y, v_z\n");

    while(ReadPositions(&read_data))
    {
        double com1_velocity[3] = {0, 0, 0}; // com velocity for polymer 1
        double com2_velocity[3] = {0, 0, 0}; // com velocity for polymer 2
        // Adding the velocities together:
        for(int i = 0; i < read_data.numberOfMonomers; i++)
        {
            if(read_data.atomTypes[i] == 1) // if atom belongs to polymer 1
            {
                for(int j = 0; j < 3; j++)
                {
                    com1_velocity[j] += read_data.bead_positions[i][j];
                }
            }
            else if(read_data.atomTypes[i] == 2) // if atom belongs to polymer 2
            {
                for(int j = 0; j < 3; j++)
                {
                    com2_velocity[j] += read_data.bead_positions[i][j];
                }
            }
            else
            {
                printf("ERROR: The atom type %i of atom %i was not recognized! Terminating.\n", read_data.atomTypes[i], i+1);
                exit(1);
            }
        }
        // Normalizing:
        double com_velocity[3]; // relative com velocity
        for(int j = 0; j < 3; j++)
        {
            com1_velocity[j] = com1_velocity[j] / numberOfMonomers;
            com2_velocity[j] = com2_velocity[j] / numberOfMonomers;
            com_velocity[j] = com1_velocity[j] - com2_velocity[j];
        }
        
        // Printing to file:
        fprintf(writeFile, "%i, %.6lf, %.6lf, %.6lf\n", read_data.timeStep, com_velocity[0], com_velocity[1], com_velocity[2]);
    }
    fclose(writeFile);
    FreeAllocatedMemory(&read_data);
}

int main(int argc, char** argv)
{
    SetConstants(argc, argv);
    CalculateAndPrintVelocities(regionPair);
    // ReadAndPrintVelocities();
    return 0;
}