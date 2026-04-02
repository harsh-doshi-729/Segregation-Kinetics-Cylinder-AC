// This is a script to read the position dump file generated in a simulation and extract the positions of the monomers for a particular timestep.
// These positions are then used as the positions to create an initial_configuration file. The bonds nad velocities are taken from another initial config file.
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "../../../../../Scripts/System_File_Paths/system_file_paths.h"
#include GENERAL_FILE_METHODS

//global variables:
char* positionDumpFilePath;
char* positionDumpDirectory;
char* initialConfigTemplateFilePath; // the file path to the initial config file from which the header, bonds section, and velocity section will be read
char* destinationDirectory; // the path to the directory where the output initial config file with be written
char* writeFilePath;
char* architecture;

const int headerLength = 9; // the number of lines in the header section of a positions dump file
int positionDataInterval; // the positions are dumped every this amount of steps during the simulation
int timeStepOfInterest; // the time step of the simulation at which the positions are desired
const int numberOfPolymers = 2;
int numberOfMonomers;
int runNumber;
const int numberOfMandatoryArguments = 2;

// bool updateBoxDimensions; // a flag to indicate whether to keep the same box dimensiosn as given in the initial configuration template file or to replace it with new dimensions defined below
// double newDimensions[3][2]; // an array that contains the dimensions of the box along the three directions: [[xlo, xhi], [ylo, yhi], [zlo, zhi]]


//Sets the values for the  global variables:
void SetConstants(int argc, char** argv)
{
    // printf("Args: %i\nThreshold: %i\n", argc, numberOfMandatoryArguments + 1);
    if(argc < numberOfMandatoryArguments + 1) // number of arguments passed is less than expected
    {
        printf("Not enough arguments passed! Please pass the architecture(string),the number of monomers(int), and the run number (int) as arguments from the command line\nFor example:\n./createInitialStateExecuatable.out Arc2 200 16\n");
        exit(1);
    }
    else
    {
        //assuming atleast two arguments were passed; aprat from the default program name
        architecture = argv[1]; // the first argument after the program name
        numberOfMonomers = atoi(argv[2]); // the second argument after the program name
        //atoi(char*) is a function that takes an argument and returns the int expressed as the chaarcter string. If the characters are not a number, it returns 0
        if(numberOfMonomers == 0)
        {
            printf("Second argument \"%s\" passed was not a valid value for number of monomers! Please pass a valid integer\n", argv[2]);
            exit(1);
        }
        runNumber = GetRunIndex(argc, argv, numberOfMandatoryArguments);
    }

    SetDirectoryPath(&positionDumpDirectory, NEW_SEGREGATION, numberOfMonomers);
    SetDirectoryPath(&destinationDirectory, CHECK_CONCATENATION, numberOfMonomers);
    Append(&initialConfigTemplateFilePath, CHECK_CONCATENATION, "b200/Arc0/initial_configuration.txt", "");
    timeStepOfInterest = 10000;
    positionDataInterval = 10000;
    
}

// //Sets the name (along with the prefix) for the position dump file
// void SetPositionDumpFileName(char** fileNamePointer)
// {
//     char* name;
//     if(numberOfMonomers == 200)
//         name = "positions.dump";
//     else if(numberOfMonomers == 500)
//         name = "visual.dump";
    
//     int bytes = asprintf(fileNamePointer, "b%i/%s/run%i/%s.dump", numberOfMonomers, architecture, runNumber, name);
//     if(bytes == -1)
//     {
//         printf("Memory could not be allocated for the positins dump file name!\n");
//         exit(1);
//     }
// }

//Sets the position file path based on the global variables
void SetReadFilePath(void)
{
    char* fileName;
    if(numberOfMonomers == 200)
        fileName = "positions.dump";
    else if(numberOfMonomers == 500)
        fileName = "visual.dump";

    char* folderName;

    SetFolderName(&folderName, architecture, runNumber);
    Append(&positionDumpFilePath, positionDumpDirectory, folderName, fileName);
    free(folderName);
}

//Sets the initial config template file path
void SetInitialTemplateFilePath(void)
{
    char* fileName = "initial_configuration.txt";
    char* folderName;
    SetFolderName(&folderName, architecture, runNumber);
    Append(&initialConfigTemplateFilePath, positionDumpDirectory, folderName, fileName);
    free(folderName);
}

//sets the filepath for the output initial configuration file
void SetWriteFilePath(void)
{
    char* fileName = "initial_configuration.txt";
    char* folderName;
    SetFolderName(&folderName, architecture, runNumber);
    Append(&writeFilePath, destinationDirectory, folderName, fileName);
    free(folderName);
}

//Frees the allocated memory for the file paths
void FreeFilePaths(void)
{
    free(positionDumpFilePath);
    free(initialConfigTemplateFilePath);
    free(writeFilePath);
}

//Copies the header section from the initial template file to the output file
void CopyHeader(FILE* readFilePointer, FILE* writeFilePointer)
{
    int buffer = 200;
    char line[buffer];
    while(true)
    {
        fgets(line, buffer, readFilePointer);

        //extracting first 5 chars from line to check if it is "Atoms":
        char word[6];
        strncpy(word, line, 5);
        word[5] = '\0'; //terminating character

        //checking if the word matches:
        if(strcmp(word, "Atoms") == 0)
        {
            // Atoms section reached and title not printed
            break;
        }
        fprintf(writeFilePointer, "%s", line); // writing the read line

    } 
    //Atoms section reached
}

//Reads the atomic positions from the position dump file and writes the Atoms section for the output file
void WriteAtomsSection(FILE* readFilePointer, FILE* writeFilePointer, int timeStepOfInterest)
{
    int buffer = 200;
    char line[buffer];

    // writing title:
    fprintf(writeFilePointer, "Atoms # angle\n");
    // printing empty line before the position data
    fprintf(writeFilePointer, "\n");

    // reaching the correct line number in the position dump file:
    int lineNumber = (timeStepOfInterest / positionDataInterval - 1) * (headerLength + numberOfPolymers * numberOfMonomers) + headerLength;
    // skipping the above lines:
    for(int i = 0; i < lineNumber; i++)
        fgets(line, buffer, readFilePointer);

    // now starting to parse the position data:
    //some constants:
    int molID = 0;
    int numberOfPositionEntries = numberOfPolymers * numberOfMonomers;
    for(int i = 0; i < numberOfPositionEntries; i++)
    {
        int atomID, atomType;
        double x, y, z;
        fscanf(readFilePointer, "%i %i %lf %lf %lf\n", &atomID, &atomType, &x, &y, &z);
        //writing to the ouput file, and inserting the molID:
        fprintf(writeFilePointer, "%i %i %i %.12lf %.12lf %.12lf\n", atomID, molID, atomType, x, y, z);
    }
    // printing empty line:
    fprintf(writeFilePointer, "\n");
    //debugging:
    // printf("Last line read in the atoms section: %s\n", line);
}

//Skips a number of lines equal to the length of a section of the file 
void SkipSection(FILE* readFilePointer, bool isAtomsSection)
{
    int lengthOfSection = 1 + 2 + numberOfPolymers * numberOfMonomers; // 1 title + 2 blank lines + data for each monomer
    if(isAtomsSection)
        lengthOfSection--; // The title for the atoms section has already been read by the COpyHeader function; thus removing one line

    int buffer = 200;
    char line[buffer];
    for(int i = 0; i < lengthOfSection; i++)
    {
        fgets(line, buffer, readFilePointer);
        //debugging:
        // if(i == 0) // first line
            // printf("First line skipped: %s\n", line);
    }
    // printf("Last line read while skipping section: %s\n", line);
}   

//Copies all the lines below the read file pointer and appends them to the output file
void CopyRestOfFile(FILE* readFilePointer, FILE* writeFilePointer)
{
    int buffer = 200;
    char line[buffer];

    while(true)
    {
        fgets(line, buffer, readFilePointer); // reading the line
        // if the end of file is reached and there is no file to read, the feof method would return a non zero value now
        //checking for End of File:
        if(feof(readFilePointer) != 0)
            break;

        //else, print the line to the writeFile
        fprintf(writeFilePointer, "%s", line);
    }
    //EOF reached
}

// Creates the new initial configuration file based on the positions at a particular timeStepand the template initial file
void CreateNewInitialConfigurationFile(int timeStepOfInterest, bool skipVelocitiesSection)
{
    // opening files:
    FILE* positionDumpFile;
    FILE* templateInitialFile;
    FILE* outputInitialFile;

    if((positionDumpFile = fopen(positionDumpFilePath, "r")) == NULL)
    {
        printf("The positions dump file could not be read!\n Path: %s\n", positionDumpFilePath);
        exit(1);
    }
    if((templateInitialFile = fopen(initialConfigTemplateFilePath, "r")) == NULL)
    {
        printf("The initial configuration template file could not be read!\n Path: %s\n", initialConfigTemplateFilePath);
        exit(1);
    }
    if((outputInitialFile = fopen(writeFilePath, "w")) == NULL)
    {
        printf("The output initial configuration file could not be opened!\n Path: %s\n", writeFilePath);
        exit(1);
    }

    CopyHeader(templateInitialFile, outputInitialFile);
    WriteAtomsSection(positionDumpFile, outputInitialFile, timeStepOfInterest);
    SkipSection(templateInitialFile, true); // skipping the atoms section
    if(skipVelocitiesSection)
        SkipSection(templateInitialFile, false); // skipping velocities section
    CopyRestOfFile(templateInitialFile, outputInitialFile);

    fclose(positionDumpFile);
    fclose(templateInitialFile);
    fclose(outputInitialFile);

    FreeFilePaths();

    printf("The initial configuration file was successfully copied and modified!\n");
}

int main(int argc, char** argv)
{
    SetConstants(argc, argv);
    // printf("Run: %i\n", runNumber);
    if(runNumber == -1)
    {
        printf("No run number was given while invoking the script. Terminating\n");
        exit(1);
    }
    SetReadFilePath();
    SetInitialTemplateFilePath();
    SetWriteFilePath();

    bool skipVelocitiesSection = false;
    CreateNewInitialConfigurationFile(timeStepOfInterest, skipVelocitiesSection);
    return 0;
}