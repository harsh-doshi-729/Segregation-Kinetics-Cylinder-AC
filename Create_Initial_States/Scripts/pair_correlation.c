// This script is used to calculate the pair correlation function for monomers part of polymers in a cylinder.
// The pair correlation function can be specified to be of all pairs, or pairs of same type, or pairs of opposite type

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <time.h>

// importing paths file:
#include "../../../../../Scripts/System_File_Paths/system_file_paths.h"
#include GENERAL_FILE_METHODS
#include LAMMPS_POSITION_FILE
#include MONOMER_DISTRIBUTION
#include CROSS_LINKS_DATABASE
#include PAIR_CORRELATION

// Global Variables:

// System information:
int numberOfMonomers; // The number of monomers in a single polymer
int numberOfPolymers = 2; // The number of polymers in the system
int totalMonomers; // The total number of monomers in the system
char* architecture; // The name of the architecture

// Simulation information:
int runIndex;
double boxVolume; // The total cylindrical volume of the box
bool isSystemGas = false; // A flag to indicate whether the system is a gas of monomers/particles in a cubical box
// This system can be used for verifying the validity of the g(r) calculation in this script
double cubeLength = 20; // The length of the cube in which the gas is contained, in case the system is a gas

// Reading file tools:
char* directory; // The path to the folder where the simulation data is stored
simulation_data read_data; // The struct to read a LAMMPS Position fump file

// Monomer distribution objects:
bool calculatePairCorrelationSeparately; // A flag to indicate whether the pair correlation function should be separately calculated
// for monomers of the same polymer and monomer pairs from different polymers along with the total; if false, it calculates g(r) only for the entire system
monomer_distribution_data all_dist_data; // The distance distribution of all pairs of monomers
monomer_distribution_data same_dist_data; // The distance distribution of pairs belonging to the same polymer
monomer_distribution_data diff_dist_data; // The distance distribution of pairs belonging to different polymers
double cutoff; // The cutoff distance beyond which the distribution should not be calculated
double binWidth = 0.2; // Default binwidth
// TODO: Add distance distributions for cross pair correlation, and same type pair correlation

// Directory:
char* directory;
char* acceptedDirectoryLabels[2] = {"Create_Initial_States", "new_segregation"};
int labelArrayLength = 2;
bool forInitialization = false; // A flag to indicate whether the the pair correlation is to be calculated for a segregation run

// Ignoring bonded monomers:
bool ignoreBondedNeighbours = true;  // A flag to indicate whether the bonded neighbours of a monomer should be ignored from the list of neighbours while calculating the pair_correlation 
// TODO: Accept the value for this flag from the user
Architecture archCrossLinks; // object to store cross links data

// Calculates and returns the volume of the cube containing the gas of monomers
double CalculateGasBoxVolume(void)
{
    return pow(cubeLength, 3);
}

// Calculates and returns the volume of the box based on the number of monomers
double CalculateBoxVolume(int numberOfMonomers)
{
    // Reading diameter from database:
    InitializeArcAndReadDiameter(&archCrossLinks, numberOfMonomers, architecture);
    double radius = archCrossLinks.confinementDiameter / 2; // radius of the cylinder
    FreeArchitecture(&archCrossLinks);
    double length = 10 * radius;
    return M_PI * pow(radius, 2) * length;
}

// Reads the passed argument as a boolean value and returns it
bool ReadFlag(char* args)
{
    if(strcasecmp(args, "false") == 0 || strcmp(args, "0") == 0)
            return false;
    else if(strcasecmp(args, "true") == 0 || strcmp(args, "1") == 0)
        return true;
    else
    {
        printf("The argument for a flag (%s) was not recognized as bool value! Terminating.\n", args);
        exit(1);
    }
}

// Sets some of the global variables with the values passed to the command line
void SetConstants(int argc, char** argv)
{
    int numberOfMandatoryArguments = 7;


    if(argc < numberOfMandatoryArguments + 1)
    {
        printf("Not enough arguments passed! Please pass the number of monomers, the name of the architecture, the run index, the cutoff distance, the bin width, the flag indicating if pair correlations should be calculated separately, and the flag indicating whether bonded monomer pairs should be ignored.\n");
        printf("If the first flag is 1 (true), pair correlation for monomers of the same polymer and pair correlation for monomers from different poymers will be calculated separately along with the total pair correlation.\n");
        printf("If the first flag is 0 (false), the pair correlation only for all monomers in the simulation will be calculated.\n");
        printf("If the second flag is 1 (true), the same and diff pair coorelations will be calculated without considering pairs of bonded monomers (cross-links, bonds within a polymer).\n");
        printf("If the second flag is 0 (false), the pair correlation is calculated in its usual way.\n");
        printf("Format: ./a.out <numberOfMonomers> <architecture> <runIndex> <cutoff> <binWidth> <pair-flag> <bond-flag>\n");
        printf("An optional argument for the calculation directory label can be passed to specify where the pair correlation is to be calculated.\n");
        printf("Accepted directory labels: ");
        PrintArray(labelArrayLength, acceptedDirectoryLabels);
        // TODO: binwidth?
        // TODO: Add optional argument to decide which pair correlation to calculate.
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

        if(isSystemGas)
            totalMonomers = 1200;
        else
            totalMonomers = numberOfPolymers * numberOfMonomers;

        cutoff = atof(argv[4]);
        if(cutoff == 0)
        {
            printf("The argument passed for cutoff distance %s could not be converted to a float! Terminating.\n", argv[4]);
            exit(1);
        }

        binWidth = atof(argv[5]);
        if(binWidth == 0)
        {
            printf("The argument passed for bin width %s could not be converted to a float! Terminating.\n", argv[5]);
            exit(1);
        }

        calculatePairCorrelationSeparately = ReadFlag(argv[6]);
        ignoreBondedNeighbours = ReadFlag(argv[7]);

        if(isSystemGas)
            boxVolume = CalculateGasBoxVolume();
        else
            boxVolume = CalculateBoxVolume(numberOfMonomers);

        // directory:
        // Checking for an optional argument:
        if(argc > numberOfMandatoryArguments + 1)
        {
            char* label = argv[numberOfMandatoryArguments+1];
            int index = SearchStringArray(label, labelArrayLength, acceptedDirectoryLabels);
            if(index == -1)
            {
                printf("The directory label %s was not recongnized. Please try again with one of the accepted labels: ", label);
                PrintArray(labelArrayLength, acceptedDirectoryLabels);
                exit(1);
            }
            char* filePrefix;
            if(index == 0) // First element: Create_Initial_States passed
                filePrefix = CREATE_INITIAL_STATES;
            else if(index == 1) // Second element: new_segregation passed
            {
                filePrefix = NEW_SEGREGATION;
                forInitialization = true;
            }

            SetDirectoryPath(&directory, filePrefix, numberOfMonomers);
        }
        else // No optional arguments passed
            SetDirectoryPath(&directory, CREATE_INITIAL_STATES, numberOfMonomers); // default directory
    }
}

// Sets the file path to the position dump file
void SetReadFilePath(char** filePathPointer)
{
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    char* fileName;
    if(isSystemGas)
        fileName = "cube_visual.dump";
    else if(forInitialization)
        fileName = "visual.dump";
    else
        fileName = "distribution_positions.dump";
    Append(filePathPointer, directory, folder, fileName);
    free(folder);
}

// Sets the complete fileName in this format: "<folderName>/<fileName>_b<binwidth><ext>"; only the fileName needs to be passed
void SetCompleteFileName(char** pointer, char* fileName)
{
    char* folderName = "pair_correlation/";
    char* extension = ".csv";
    int bytes = asprintf(pointer, "%s%s_b%02i%s", folderName, fileName, (int)round(binWidth * 10), extension);
    if(bytes == -1)
    {
        printf("Memory could not be allocated for the file name! Terminating.\n");
        exit(1);
    }
}



// Sets the file path to a file where some data will be printed; the filename should reflect what the data is
void SetWriteFilePath(char** filePathPointer, char* fileName)
{
    // TODO: Accept additional arguments to specify which pair correlation is to be printed?
    char* folder;
    SetFolderName(&folder, architecture, runIndex);
    char* completeFileName;
    SetCompleteFileName(&completeFileName, fileName);
    Append(filePathPointer, directory, folder, completeFileName);
    free(completeFileName);
    free(folder);
}

// Computes and returns the distance calculated between a pair of monomers
double CalculateDistance(int index1, int index2)
{
    // Verifying if the indices are within bounds:
    int totalMonomers = read_data.numberOfMonomers;
    if(index1 >= 0 && index2 >=0 && index1 < totalMonomers && index2 < totalMonomers) // indices lie between 0 and numberOfMonomers - 1
    {
        double distance = 0;
        for(int i = 0; i < 3; i++)
        {
            // enforcing minimum image convention if the system is a gas
            double delta = fabs(read_data.bead_positions[index1][i] - read_data.bead_positions[index2][i]);
            if(isSystemGas)
            {
                if(delta > cubeLength / 2)
                    delta = cubeLength - delta;
            }
            distance += pow(delta, 2);
        }
        distance = sqrt(distance);
        return distance;
    }
    else
    {
        printf("The index %i or %i is out of bounds!\n", index1, index2);
        exit(1);
    }
}

// Determines if the passed pair of monomer indices is valid according to the flags passed;
// See docstring of AddDistancesToDIstribution() for more information on the flags
bool IsIndicesPairValid(int index1, int index2, bool fromSamePolymer, bool fromDifferentPolymers, bool ignoreBondedNeighbours)
{
    bool areIndicesFromSamePolymer = (read_data.atomTypes[index1] == read_data.atomTypes[index2]);
    // Note: the above flag is independent of the argument fromSamePolymer
    // The argument fromSamePolymer indicates the requirement of whether only indices from same polymers are valid

    // Default case of both flags true:
    if(fromSamePolymer && fromDifferentPolymers)
        return true;
    // Handling the invalid cases:
    if(fromSamePolymer && !areIndicesFromSamePolymer)
    {
        return false;
    }        
    if(fromDifferentPolymers && areIndicesFromSamePolymer)
    {
        return false;
    }
    // Ignoring pair of bonded monomers if the ignoreBondedNeighbours flag is true
    if(ignoreBondedNeighbours)
    {
        if(AreMonomersBonded(archCrossLinks, numberOfMonomers, index1, index2))
            return false; // ignoring such bonded monomers
    }
    // else: all other cases are valid
    return true;
}


// Adds the distances between pairs of monomers calculated from a single timestep to the distance distribution
// Takes the distribution_data struct, and two flag. The first flag indicates whether only pairs of monomers from the same polymer should be considered.
// The second flag indicates whether only pairs of monomers from different polymers should be considered. The flags are mututally exclusive
// The third flag indicates whether the distances between bonded monomers should not be added for the pair_correlation calculation
void AddDistancesToDistribution(monomer_distribution_data* dist_data, bool fromSamePolymer, bool fromDifferentPolymers, bool ignoreBondedNeighbours)
{
    if(!fromSamePolymer && !fromDifferentPolymers)
    {
        printf("AddDistanceToDistribution(): Both flags cannot be false! Terminating.\n");
        exit(1);
    }
    // For all pairs of monomers; regardless of type
    for(int i = 0; i < read_data.numberOfMonomers; i++)
    {
        for(int j = i + 1; j < read_data.numberOfMonomers; j++)
        {
            if(IsIndicesPairValid(i, j, fromSamePolymer, fromDifferentPolymers, ignoreBondedNeighbours))
            {    
                double distance = CalculateDistance(i, j);
                if(distance > cutoff)
                    continue;
                IncrementMonomerCounts(dist_data, distance);
                IncrementMonomerCounts(dist_data, distance); // Here, the distance is the position coordinate of the distribution;
                // Each count is not a monomer, but a pair of monomers
                // Adding each distance twice since the distance is symmetric with respect to interchange of monomer indices;
                // both ij and ji distances taken into account for normalization
            }
        }
    }
}

// Prints the distance distribution to a file
void PrintDistanceDistribution(FILE* writeFilePointer, monomer_distribution_data dist_data)
{
    // printing header:
    fprintf(writeFilePointer, "Left Bin Edge, Probability distribution for distance\n");
    
    // Printing the distance distribution after proper normalization:
    double normalizationFactor = dist_data.numberOfDataPoints * dist_data.binWidth;
    // I am not using the Monomer_Distribution method to normalize a distribution since I might need the unnormalized distribution later
    for(int i = 0 ; i < dist_data.numberOfBins; i++)
    {
        double leftBinEdge = i * dist_data.binWidth;
        fprintf(writeFilePointer, "%lf, %lf\n", leftBinEdge, dist_data.monomerCounts[i] / normalizationFactor);
    }
}

// Prints the pair correlation function by properly normalizing the distance distribution
// Normalizing done with the average density and number of steps the data is collected
void PrintPairCorrelation(FILE* writeFilePointer, monomer_distribution_data dist_data, int numberOfTimeSteps)
{
    // printing header:
    fprintf(writeFilePointer, "Left Bin Edge, Pair Correlation Function\n");

    // g(r) = 1/N * dN(r, r+dr) / (4 PI r^2 dr <rho>)
    // And then normalizing with the number of time steps
    double averageDensity = dist_data.numberOfMonomers / boxVolume;
    // The normalization factor without the radius^2 in the denominator; common for all radii
    double preNormalization = dist_data.numberOfMonomers * 4 * M_PI * dist_data.binWidth * averageDensity * numberOfTimeSteps;
    for(int i = 0; i < dist_data.numberOfBins; i++)
    {
        double radius = i * binWidth; // The radius of the shell being considered for g(r); aka left bin edge
        double normalizationFactor = preNormalization;
        if(i != 0) // Not normalizing by zero radius;ignoring the radius term
            normalizationFactor = normalizationFactor * pow(radius, 2);
        fprintf(writeFilePointer, "%lf, %lf\n", radius, dist_data.monomerCounts[i] / normalizationFactor);
    }
}

// Calculates the distance distribution(s) for the entire simulation 
void ComputeEntireDistribution(void)
{
    // Setting up reading file:
    char* readFilePath;
    SetReadFilePath(&readFilePath); // memory will be handled by the read_data object
    InitializeSimulationData(&read_data, totalMonomers, readFilePath);

    // Setting up distribution:
    bool isBoxCentredAtOrigin = false; // The position coordinate is distance which is always positive;
    // thus we want the "box" to start from 0
    // The number of monomers passed in the initialization method is the number of other monomers in the box that will be seen by the central atom of g(r)
    // For all pair correlation functions, a common number equal to the total number of monomers in the cylinder has been taken
    // This is so that the separate pair correlations add to give the total pair correlation
    if(calculatePairCorrelationSeparately)
    {
        InitializeDistributionDataWithAxisLengthBinWidth(&same_dist_data, totalMonomers, cutoff, isBoxCentredAtOrigin, binWidth);
        InitializeDistributionDataWithAxisLengthBinWidth(&diff_dist_data, totalMonomers, cutoff, isBoxCentredAtOrigin, binWidth);
    }
    if(ignoreBondedNeighbours || !calculatePairCorrelationSeparately)  // calculating total monomers if the bonded monomers are being ignored or if not calculating separately
    {
        // total pair correlation is calculated in all cases:
        // It is is explicitly calculated when bonded monomers are being ignore even if the separate ones are being calculated
        // because the total pair correlation cannot be reconstructed by adding the two components
        InitializeDistributionDataWithAxisLengthBinWidth(&all_dist_data, totalMonomers, cutoff, isBoxCentredAtOrigin, binWidth);
    }
    // Setting up the write file for pair correlation:
    char* writeFilePaths[3];
    int numberOfFiles = 0; // The number of pair correlation functions to be written
    // default case:
    numberOfFiles = 1;
    SetWriteFilePath(&(writeFilePaths[0]), "pair_correlation");
    if(calculatePairCorrelationSeparately)
    {
        numberOfFiles = 3;
        SetWriteFilePath(&(writeFilePaths[1]), "same_pair_correlation");
        SetWriteFilePath(&(writeFilePaths[2]), "diff_pair_correlation");
    }
    FILE* pairFilePointers[3];
    for(int i = 0; i < numberOfFiles; i++)
    {
        pairFilePointers[i] = fopen(writeFilePaths[i], "w"); // Opening the write file
        if(pairFilePointers[i] == NULL)
        {
            printf("The write file %s could not be opened! Terminating.\n", writeFilePaths[i]);
            exit(1);
        }
        free(writeFilePaths[i]);
    }
    // Setting up the distance distribution write file: debugging
    // SetWriteFilePath(&writeFilePath, "distance_distribution");
    // FILE* distFilePointer = fopen(writeFilePath, "w");
    //     if(distFilePointer == NULL)
    // {
    //     printf("The write file %s could not be opened! Terminating.\n", writeFilePath);
    //     exit(1);
    // }
    // free(writeFilePath);

    // Initializing Cross Links object in case the bonded monomers need to be ignored
    InitializeArcAndReadCrossLinks(&archCrossLinks, numberOfMonomers, architecture);
    // Reading the position dump file a timestep at a time:
    int numberOfTimeSteps = 0; // Counter to count the number of microstates taken to calculate pair correlation
    while(ReadPositions(&read_data))
    {
        if(calculatePairCorrelationSeparately)
        {
            AddDistancesToDistribution(&same_dist_data, true, false, ignoreBondedNeighbours); // Same polymer pair_correlation
            AddDistancesToDistribution(&diff_dist_data, false, true, ignoreBondedNeighbours); // Different polymers pair_correlation
        }
        if(ignoreBondedNeighbours || !calculatePairCorrelationSeparately) // explicitly calculating total monomers if the bonded monomers are being ignored or if not calculating separately
            AddDistancesToDistribution(&all_dist_data, true, true, false); // All pair of monomers in the system; even bonded monomers
        numberOfTimeSteps++;
    }

    if(calculatePairCorrelationSeparately && !ignoreBondedNeighbours) // can reconstruct the total pair correaltion by summing two components
    {
        // Adding the two separate distributions to get the total:
        AddMonomerCounts(same_dist_data, diff_dist_data, &all_dist_data); // all_dist_data will be initialized
        all_dist_data.numberOfMonomers = totalMonomers; // Updating number of monomers to be the total monomers in the system
    }
    // Freeing cross links memory:
    FreeArchitecture(&archCrossLinks);

    // Print the distance distribution: debugging
    // PrintDistanceDistribution(distFilePointer, dist_data);
    // fclose(distFilePointer);
    
    // Printing the correlation function(s):
    monomer_distribution_data dist_data_objects[3] = {all_dist_data, same_dist_data, diff_dist_data};
    for(int i = 0; i < numberOfFiles; i++)
    {
        PrintPairCorrelation(pairFilePointers[i], dist_data_objects[i], numberOfTimeSteps);
        fclose(pairFilePointers[i]);
    }

    // Printing number of timesteps read to the console:
    printf("Number of time steps from which data was collected: %i\n", numberOfTimeSteps);

    // Freeing memory:
    FreeAllocatedMemory(&read_data);
    for(int i = 0 ; i < numberOfFiles; i++)
        FreeMonomerDistribution(&dist_data_objects[i]);
}

// Calculates the total pair correlation function using the new Generalised Particles definition for consistency check
void ComputeGenParticlePairCorrelation(void)
{
    // Initializing read file:
    char* readFilePath;
    SetReadFilePath(&readFilePath);
    InitializeSimulationData(&read_data, totalMonomers, readFilePath);
    free(readFilePath);

    // Initializing monomer distribution (histogram):
    monomer_distribution_data histogram;
    InitializeDistributionDataWithAxisLengthBinWidth(&histogram, totalMonomers, cutoff, false, binWidth);

    // Setting up an array of Generalised particles:
    GeneralisedParticle genParticles[totalMonomers];

    // Reading file:
    int numberOfTimeSteps = 0;
    while(ReadPositions(&read_data))
    {
        // Assigning positions to generalised particles:
        for(int i = 0; i < totalMonomers; i++)
            UpdateParticleCoordinates(&(genParticles[i]), read_data.bead_positions[i]);
        UpdateDistanceHistogram(&histogram, genParticles);
        numberOfTimeSteps++;
    }
    // calculating the pair correlation:
    double pairCorrelation[histogram.numberOfBins];
    CalculatePairCorrelationFunction(&histogram, pairCorrelation, boxVolume, numberOfTimeSteps);

    // Printing the pair correlation:
    char* writeFilePath;
    SetWriteFilePath(&writeFilePath, "pair_correlation");
    FILE* writeFilePointer = fopen(writeFilePath, "w");
    if(writeFilePointer == NULL)
    {
        printf("The file '%s' could not be opened for writing!\n", writeFilePath);
        exit(1);
    }
    free(writeFilePath);
    // Writing to file:
    fprintf(writeFilePointer, "Bin Left Edge, Pair Correlation Function\n"); // header
    for(int i = 0; i < histogram.numberOfBins; i++)
        fprintf(writeFilePointer, "%0.6lf, %0.6lf\n", i * binWidth, pairCorrelation[i]);
    fclose(writeFilePointer);

    // Freeing memory:
    FreeAllocatedMemory(&read_data);
    FreeMonomerDistribution(&histogram);
}

int main(int argc, char** argv)
{
    SetConstants(argc, argv);
    
    // Timing the calculation:
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);
    printf("Start: %s\n", asctime(tm));
    clock_t start = clock();

    ComputeEntireDistribution();
    free(directory);

    t = time(NULL);
    tm = localtime(&t);
    clock_t finish = clock();
    long int elapsedTime = (double)(finish - start)/ CLOCKS_PER_SEC;
    printf("Finish: %s\nElapsed (CPU) Time: %lih %lim %lis\n", asctime(tm), elapsedTime/3600, (elapsedTime/60) % 60, elapsedTime%60);

    return 0;
}