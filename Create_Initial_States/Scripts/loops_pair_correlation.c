// This script used the Generalised Particle pair correlation script to calculate the pari correlation function for different regions of a polymer
// The cross-links demarcate the regions, and each region(usually a loop) is abstracted as a particle with its position being the CoM of the region

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <time.h>

// importing paths file:
#include "../../../../../Scripts/System_File_Paths/system_file_paths.h"
#include GENERAL_FILE_METHODS
#include LAMMPS_POSITION_FILE
#include CROSS_LINKS_DATABASE
#include MONOMER_DISTRIBUTION
#include PAIR_CORRELATION

// Global variables:
int numberOfMonomers;
int numberOfPolymers = 2;
int totalMonomers;
char* architecture;
int runIndex = -1;

// Directory:
char* directory; // The path to the directory where the position dump file exists
int labelArrayLength = 2;
char* acceptedDirectoryLabels[2] = {"Create_Initial_States", "new_segregation"};
bool forInitialization = true; // A flag to indicate whether the operations are to be performed for the initialization simulations; default is true

// Pair Correlation function details:
double binWidth; // The binwidth for the pair correlation function
double cutoff; // The cutoff distance up to which the pair correlation function calculation will be restricted

// Reads the command line arguments and sets the required constants
void SetConstants(int argc, char** argv)
{
    int numberOfMandatoryArguments = 5;
    if(argc < numberOfMandatoryArguments + 1) // Not enough arguments passed
    {
        printf("Not enough arguments passed! Please pass the number of monomers (integer), the architecture name (string), the run index (integer), the cutoff distance (float), and the bin width (float) as arguments to the command line while running the executable in the following format:\n");
        printf("./pair.out <numberOfMonomers> <architecture> <runIndex> <cutoff> <binWidth>\n");
        printf("An optional argument for the directory label can be passed to perform the operations in that particular directory.\n");
        printf("Accepted directory labels: ");
        PrintArray(labelArrayLength, acceptedDirectoryLabels);
        exit(1);
    }
    else
    {
        // Assigning variables:
        numberOfMonomers = atoi(argv[1]);
        if(numberOfMonomers == 0)
        {
            printf("ERROR: The passed argument '%s' is not a valid number of monomers.\n", argv[1]);
            exit(1);
        }
        totalMonomers = numberOfMonomers * numberOfPolymers;

        architecture = argv[2];

        runIndex = atoi(argv[3]);
        if(runIndex == 0)
        {
            printf("ERROR: The passed argument '%s' is not a valid run index.\n", argv[3]);
            exit(1);
        }

        cutoff = atof(argv[4]);
        if(cutoff == 0)
        {
            printf("ERROR: The passed argument '%s' is not a valid cutoff distance.\n", argv[4]);
            exit(1);
        }
        binWidth = atof(argv[5]);
        if(binWidth == 0)
        {
            printf("ERROR: The passed argument '%s' is not a valid bin width.\n", argv[5]);
            exit(1);
        }

        // Optional argument(s):
        if(argc > numberOfMandatoryArguments + 1) // Optional argument(s) passed
        {
            char* directoryLabel = argv[numberOfMandatoryArguments + 1];
            char* filePrefix;
            int index = SearchStringArray(directoryLabel, labelArrayLength, acceptedDirectoryLabels);
            if(index == -1)
            {
                printf("ERROR: The optional argument '%s' was not recognized as a directory label. Please pass one of the following for the directory label: ", directoryLabel);
                PrintArray(labelArrayLength, acceptedDirectoryLabels);
                exit(1);
            }
            else if(index == 0) // Create_Initial_States passed
            {
                filePrefix = CREATE_INITIAL_STATES;
                forInitialization = true;
            }
            else if(index == 1) // new_segregation passed
            {
                filePrefix = NEW_SEGREGATION;
                forInitialization = false;
            }
            // Setting directory:
            SetDirectoryPath(&directory, filePrefix, numberOfMonomers);
        }
        else // No optional arguments passed
            SetDirectoryPath(&directory, CREATE_INITIAL_STATES, numberOfMonomers); // default directory
    }
}

// int test_main(int argc, char** argv)
// {
//     // testing the SetConstants() function and checking file paths:
//     SetConstants(argc, argv);
//     printf("Number of Monomers: %i\n", numberOfMonomers);
//     printf("Architecture: %s\n", architecture);
//     printf("Run Index: %i\n", runIndex);
//     printf("Cutoff Distance: %.4lf\n", cutoff);
//     printf("Bin Width: %.4lf\n", binWidth);
//     printf("Directory: %s\n", directory);
//     char* readFilePath;
//     SetPositionDumpFilePath(&readFilePath, directory, architecture, runIndex, forInitialization);
//     printf("Read File Path: %s\n", readFilePath);
//     free(readFilePath);
//     char* writeFilePath;
//     SetCustomWriteFilePath(&writeFilePath, directory, architecture, runIndex, "loop_pair_correlation.csv");
//     printf("Write File Path: %s\n", writeFilePath);
//     free(writeFilePath);
// }

// Initializes the Generalised Particles for the loops of the two polymers
// Args:
// - genParticles: A 2D array containing the GeneralisedParticle instances for each loop of the two polymers; 1st index: index of polymer, 2nd index: index of loop
void InitializeGenParticlesArray(Architecture archCrossLinks, GeneralisedParticle genParticles[numberOfPolymers][archCrossLinks.numberOfCrossLinks])
{
    int uid = 1;
    for(int i = 0; i < numberOfPolymers; i++)
    {
        for(int j = 0; j < archCrossLinks.numberOfCrossLinks; j++)
        {
            char* description;
            int bytes = asprintf(&description, "Polymer %i, Loop %i", i+1, j+1);
            if(bytes == -1)
            {
                printf("ERROR: Memory could not be allocated for the description of the Generalised particle 'Polymer %i, Loop %i'.\n", i+1, j+1);
                exit(1);
            }
            genParticles[i][j] = GetNewGenParticle(uid++, description);
            free(description);
        }
    }
}

// Calculates the centre of mass of a particular loop of a polymer
// Args:
// - centerOfMass: The array (vector) of the 3D center of mass where the calculated values will be stored
// - archCrossLinks: The Architecture struct that contains the information about the cross-links
// - polymerIndex: The index of the polymer of which the loop is a part of; indexed from 0
// - loopIndex: The index of the loop within the polymer; indexed from 0
// - read_data: The simulation_data struct that contains the position information of the current timestep
void CalculateLoopCenterOfMass(double centerOfMass[3], Architecture archCrossLinks, int polymerIndex, int loopIndex, simulation_data read_data)
{
    // Initializing center of mass:
    for(int i = 0; i < 3; i++)
        centerOfMass[i] = 0;
    int counter = 0; // To count the number of monomers in the loop region
    // Getting the monomer index (index from 1) ast the start of the cross-link region
    int monomerIndex = archCrossLinks.crossLinks[loopIndex][0] + polymerIndex * numberOfMonomers;
    int endIndex = GetNextMonomer(archCrossLinks.crossLinks[loopIndex][1], numberOfMonomers) + polymerIndex * numberOfMonomers; // The index of monomer such beyond the end of the loop region
    // moving through the contiguous region till the end:
    do
    {
        for(int i = 0; i < 3; i++)
            centerOfMass[i] += read_data.bead_positions[monomerIndex-1][i]; // monomerIndex is 1-based
        counter++;
        // Iterating to next monomer:
        monomerIndex = GetNextMonomer(monomerIndex, numberOfMonomers);
    } while (monomerIndex != endIndex);
    // Sanity check for number of monomers in loop region
    if(counter != GetCrossLinkRegionCount(archCrossLinks.crossLinks[loopIndex], numberOfMonomers))
    {
        printf("ERROR: The number of monomers counted in the loop %i of Polymer %i does not match the number given by Cross-Links_Database calculation.\n", loopIndex, polymerIndex);
        exit(1);
    }
    // Completing the calculation for center of mass:
    for(int i = 0; i < 3; i++)
        centerOfMass[i] /= counter;
    
}

// Prints the coordinates of the Generalised Particles to the file for the current timestep
void PrintGeneralisedParticles(FILE* filePointer, int timeStep, Architecture archCrossLinks, GeneralisedParticle genParticles[numberOfPolymers][archCrossLinks.numberOfCrossLinks])
{
    for(int i = 0; i < numberOfPolymers; i++)
    {
        for(int j = 0; j < archCrossLinks.numberOfCrossLinks; j++)
        {
            fprintf(filePointer, "%i, %i, %i, %.6lf, %.6lf, %.6lf\n", timeStep, i+1, j+1, genParticles[i][j].coordinates[0], genParticles[i][j].coordinates[1], genParticles[i][j].coordinates[2]);
        }
    }
}

// Prints the pair correlation function passed to a file
// Args:
// - numberOfBins: The number of bins over which the pair correlation function has been calculated
// - binWidth: The bin width used for the histogram and pair correlation function
// - pairCorrelation: The double array that contains the values of the pair correlation function at each distance
void WritePairCorrelationFunction(int numberOfBins, double binWidth, double pairCorrelation[numberOfBins])
{
    // Setting up the write file name:
    char* writeFileName;
    int bytes = asprintf(&writeFileName, "pair_correlation/loops_pair_correlation_b%02i.csv", (int) round(binWidth * 10));
    if(bytes == -1)
    {
        printf("ERROR: Memory could not be allocated for the the loops pair correlation write file name.\n");
        exit(1);
    }
    char* writeFilePath;
    SetCustomWriteFilePath(&writeFilePath, directory, architecture, runIndex, writeFileName);
    free(writeFileName);
    FILE* writeFilePointer = fopen(writeFilePath, "w");
    if(writeFilePointer == NULL)
    {
        printf("ERROR: The file '%s' could not be opened to write the pair correlation function.\n", writeFilePath);
        exit(1);
    }
    free(writeFilePath);
    
    // Writing header:
    fprintf(writeFilePointer, "Bin Left Edge, Pair Correlation Function\n");
    // Writing values:
    for(int i = 0; i < numberOfBins; i++)
        fprintf(writeFilePointer, "%.6lf, %.6lf\n", i * binWidth, pairCorrelation[i]);
    fclose(writeFilePointer);
}

// Calculates the pair correlation function for the regions (loops) demarcated by the cross-links of the architecture.
// Writes the computed pair correlation function to a file
void ComputeLoopsPairCorrelation(void)
{
    // Setting up reading file machinery:
    char* readFilePath;
    SetPositionDumpFilePath(&readFilePath, directory, architecture, runIndex, forInitialization);
    simulation_data read_data;
    InitializeSimulationData(&read_data, totalMonomers, readFilePath);
    free(readFilePath);

    // Reading cross-links for the architecture:
    Architecture archCrossLinks;
    InitializeArcAndReadCrossLinks(&archCrossLinks, numberOfMonomers, architecture);

    // Setting up histogram:
    bool isCentredAtOrigin = false; // The pair correlation function only deals with positive distances; all values in histogram are above 0
    int numberOfGenParticles = archCrossLinks.numberOfCrossLinks; // The number of Generalised Particles (loops) for a single polymer
    monomer_distribution_data histogram; // The data from both polymers will be added to this histogram
    InitializeDistributionDataWithAxisLengthBinWidth(&histogram, numberOfGenParticles, cutoff, isCentredAtOrigin, binWidth);

    // Initializing the General Particles for the loops of the two polymers:
    GeneralisedParticle genParticles[numberOfPolymers][archCrossLinks.numberOfCrossLinks];
    InitializeGenParticlesArray(archCrossLinks, genParticles);

    // File to write the positions of Generalised particles for a sanity check:
    // char* genParticleFilePath;
    // SetCustomWriteFilePath(&genParticleFilePath, directory, architecture, runIndex, "pair_correlation/loop_com.csv");
    // FILE* genParticleFile = fopen(genParticleFilePath, "w");
    // if(genParticleFile == NULL)
    // {
    //     printf("ERROR: The file '%s' could not be opened for writing.\n", genParticleFilePath);
    //     exit(1);
    // }
    // free(genParticleFilePath);
    // // Printing header:
    // fprintf(genParticleFile, "TimeStep, Polymer Index, Loop Index, CoM_x, CoM_y, CoM_z\n");


    // Reading positions and updating histogram:
    int numberOfDataSets = 0; // The number of sets of position data, each consisting of the positions of the loops of one polymer; the ffective number of timesteps for histogram
    while(ReadPositions(&read_data))
    {
        // Updating generalised particle positions; for each loop:
        for(int i = 0; i < numberOfPolymers; i++)
        {
            for(int j = 0; j < archCrossLinks.numberOfCrossLinks; j++)
            {
                double centerOfMass[3];
                CalculateLoopCenterOfMass(centerOfMass, archCrossLinks, i, j, read_data);
                UpdateParticleCoordinates(&(genParticles[i][j]), centerOfMass);
            }
        }
        // Printing the centers of mass of each loop to a file for a sanity check: to verify if the calculation is correct
        // if(read_data.timeStep % 100000 == 0)
        //     PrintGeneralisedParticles(genParticleFile, read_data.timeStep, archCrossLinks, genParticles);

        // Adding to the pair correlation distribution below:
        for(int i = 0; i < numberOfPolymers; i++)
        {
            UpdateDistanceHistogram(&histogram, genParticles[i]); // Adding the distances between the loop gen particles for the ith polymer
            numberOfDataSets++;
        }
    }
    // fclose(genParticleFile);
    // Calculating the pair correlation function by normalising the histogram:
    double pairCorrelation[histogram.numberOfBins];
    // Reading diameter and axis length of confining cylinder; to set the box volume:
    ReadDiameterAndAxisLength(&archCrossLinks, numberOfMonomers);
    double boxVolume = CalculateConfinementVolume(archCrossLinks);
    CalculatePairCorrelationFunction(&histogram, pairCorrelation, boxVolume, numberOfDataSets); // The numberOfDataSets argument has been passed as the numberOfTimeSteps 
    // Writing the pair correlation function to a file:
    WritePairCorrelationFunction(histogram.numberOfBins, histogram.binWidth, pairCorrelation);

    // Freeing memory:
    FreeAllocatedMemory(&read_data);
    FreeMonomerDistribution(&histogram);
    for(int i = 0; i < numberOfPolymers; i++)
        for(int j = 0; j < archCrossLinks.numberOfCrossLinks; j++)
            FreeGenParticle(&(genParticles[i][j]));
    FreeArchitecture(&archCrossLinks);
}

int main(int argc, char** argv)
{
    SetConstants(argc, argv);

    // Timing the calculation:
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);
    printf("Start: %s\n", asctime(tm));
    clock_t start = clock();

    ComputeLoopsPairCorrelation();

    t = time(NULL);
    tm = localtime(&t);
    clock_t finish = clock();
    long int elapsedTime = (double)(finish - start)/ CLOCKS_PER_SEC;
    printf("Finish: %s\nElapsed (CPU) Time: %lih %lim %lis\n", asctime(tm), elapsedTime/3600, (elapsedTime/60) % 60, elapsedTime%60);
    return 0;
}
