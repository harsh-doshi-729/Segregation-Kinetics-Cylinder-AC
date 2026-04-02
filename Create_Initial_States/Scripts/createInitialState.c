// This is a script to read an initial data file for two Arc0 polymers and copy and edit it to 
// make it a valid data file for running segregation of a particular architecture.
// This is done by introducing the additional cross links for that particular architecture

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "../../Global_Scripts/System_File_Paths/system_file_paths.h"
// NOTE: Please ensure to add the option -I <path/to/project/directory> while compiling any of the C scripts that import this header file
// This makes it so that the compiler looks for the header files in the project directory as well, and thus, can find the header file paths defined in that file
#include GENERAL_FILE_METHODS
#include CROSS_LINKS_DATABASE
#include MODULO
// Config file:
#include STATE_CONFIG

// global variables:
char* directory;
char* architecture;
int numberOfPolymers;
int numberOfMonomers;
int runIndex;
int numberOfMandatoryArguments = 2;
bool useDifferentBondTypeForCrossLinks = false;
double* offsetAngles; // An array of offset angles (in degrees) read from the top of the initial file; one for each polymer

Architecture archCrossLinks; // The cross links architecture object

// Sets the values of the global constant variables based on the arguments passed
void SetConstants(int argc, char** argv)
{
    char* templateDirectory = HOME_DIR CREATE_INITIAL_STATES; // default value for the location of the initial configuration template
    if(argc < numberOfMandatoryArguments + 1) // number of arguments passed is less than expected
    {
        printf("Not enough arguments passed! Please pass the architecture(string) and the number of monomers(int) as arguments from the command line\nFor example:\n./createInitialStateExecuatable.out 200 Arc-1-2\n");
        printf("Optional arguments for the particular run number and/or a template location directory can also be passed after the mandatory arguments if desired\n");
        printf("If a run index (integer) is passed, then the configuration write file is saved to a run folder. If a template location directory is passed, the entire process is performed in that directory.\n");
        exit(1);
    }
    else
    {
        //assuming atleast two arguments were passed; aprat from the default program name
        architecture = argv[2]; // the first argument after the program name

        numberOfMonomers = atoi(argv[1]); // the second argument after the program name
        //atoi(char*) is a function that takes an argument and returns the int expressed as the character string. If the characters are not a number, it returns 0
        if(numberOfMonomers == 0)
        {
            printf("Second argument \"%s\" passed was not a valid value for number of monomers! Please pass a valid integer\n", argv[1]);
            exit(1);
        }
        
        if(argc > numberOfMandatoryArguments + 1)
        {
            runIndex = atoi(argv[numberOfMandatoryArguments+1]);
            if(argc == numberOfMandatoryArguments + 1 + 1) // exactly one extra argument
            {
                if(runIndex == 0) // validating run index
                {
                    printf("The optional argument %s could not be converted to a valid run index! Assuming it is the directory.\n", argv[numberOfMandatoryArguments+1]);
                    templateDirectory = argv[numberOfMandatoryArguments + 1];
                    // TODO: check if the directory exists
                }
                else
                {
                    printf("No argument for the optional directory passed. Continuing with default Create_Initial_States/\n");
                }
            }
            else // more than 1 argument passed
            {
                templateDirectory = argv[numberOfMandatoryArguments + 2];
                // TODO: check if the directory exists
            }
        }
    }
    // Set directory adds the corresponding subdirectory for the number of monomers to the passed folderPath
    SetDirectoryPath(&directory, templateDirectory, numberOfMonomers);
}

// Initializes the offset angles array
// If the angles are found, then it tokenizes the angles string using the delims passed
void InitializeOffsetAngles(bool areAnglesFound, char* anglesString, char* delims)
{
    offsetAngles = calloc(numberOfPolymers, sizeof(*offsetAngles));
    if(areAnglesFound)
    {
        char* token; // String to store the token extracted from the angles string
        int counter = 0;
        while((token = strtok_r(anglesString, delims, &anglesString)))
        {
            // the token should have been stored
            offsetAngles[counter++] = atof(token);
        }   
    }
}

//Sets the file name of initial configuration to be read
void SetInitialStateFileName(char** fileNamePointer)
{   
    int bytes = 0;
    if(archCrossLinks.isArcLinear)
        bytes = asprintf(fileNamePointer, "linear_initial_configuration.txt");
    else
        bytes = asprintf(fileNamePointer, "initial_configuration.txt");

    if(bytes == -1)
    {
        printf("Memory could not be allocated for the template initial file name!\n");
        exit(1);
    }
}

//Sets the read file path of the template initial configuration file
void SetReadFilePath(char** filePathPointer)
{
    char* fileName;
    SetInitialStateFileName(&fileName);
    if(USE_COMMON_TEMPLATE)
    {
        int bytes = asprintf(filePathPointer, "%sb%i/%s", HOME_DIR CREATE_INITIAL_STATES, numberOfMonomers, fileName);
        if(bytes == -1)
        {
            printf("Memory could not be allocated for the initial configuration template file name!\n");
            exit(1);
        }
    }
    else
        Append(filePathPointer, "", directory, fileName);
    free(fileName);
}

//Sets the write file path in the destination directory
void SetWriteFilePath(char** filePathPointer)
{
    char* folderName;
    char* fileName = "initial_configuration.txt";
    if(runIndex == 0) // runIndex not set
        Append(&folderName, architecture, "/", "");
    else
        SetFolderName(&folderName, architecture, runIndex);
    Append(filePathPointer, directory, folderName, fileName);
    free(folderName);
}

// Copies lines from the readFile to the writeFile till the query word is reached
void CopyLinesTillQuery(FILE* readFilePointer, FILE* writeFilePointer, char* query, bool doNotWrite)
{
    int wordLength = strlen(query);
    char word[wordLength + 1]; // looking for the word that matches the query
    int buffer = 200;
    char line[buffer];
    
    while(true)
    {
        fgets(line, buffer, readFilePointer);
        //looking for the query word in the line:
        strncpy(word, line, wordLength);
        word[wordLength] = '\0'; //terminating character
        if(strcasecmp(word, query) == 0)
            break; // does not copy the line where the query is found

        //else:
        if(!doNotWrite)
            fprintf(writeFilePointer, "%s", line); // writes the read line only if doNotWrite is false
    }
}

// Compares the query and the start of the passed string; returns true on a match and flase otherwise
bool ContainsQuery(const char* line, const char* query)
{
    int queryLength = strlen(query);
    char copy[queryLength+1]; // string to copy part of the line string
    strncpy(copy, line, queryLength); // copying the first part of 'line' to 'copy'
    // Adding the null terminating character:
    copy[queryLength] = 0;
    if(strcasecmp(copy, query) == 0)
        return true;
    else
        return false;
}

// Reads the passed string and determines if offset angles have been passed
// If the offset angle have been written, then the function reads the offset angles as strings 
// and stores them in the passed storage location
bool ReadOffsetAngles(char* line, char* offsetAnglesStorage)
{
    int lineSize = strlen(line);
    // Checking if the line contains offset angles
    char* query = "# Offset angles: ";
    if(ContainsQuery(line, query))
    {
        int queryLength = strlen(query);
        strncpy(offsetAnglesStorage, &line[queryLength], lineSize - queryLength); // copying the offset angles line after the query;
        // &line[queryLength] points to the character just after the query string in line
        offsetAnglesStorage[lineSize - queryLength - 1] = 0; // Setting the last character as the null terminating character; excluding \n
        // // Debugging information:
        // printf("Found offset angles: '%s'\n", offsetAnglesStorage);
        return true;
    }
    else
    {
        return false;
    }
    
}

//Copies the header section of the data file; correctly substitutes the bond information
void CopyHeader(FILE* readFilePointer, FILE* writeFilePointer)
{
    //copying the first few lines as is:
    // The first line is a comment explaining the config file
    // The second line is a comment containing the offset angles (optional)
    // The next line is a blank line
    int buffer = 200;
    char line[buffer];
    char anglesString[buffer];
    bool areAnglesFound = false;
    while(true)
    {
        fgets(line, buffer, readFilePointer); // reading line
        fprintf(writeFilePointer, "%s", line); // writing line

        if(!areAnglesFound)
            areAnglesFound = ReadOffsetAngles(line, anglesString);

        if(line[0] == '\n') // empty line
            break;
    }
    if(!areAnglesFound) // printing warning if angles are not found
        printf("WARNING: No offset angles found at the beginning of the file. Assuming all offset angles are 0\n");

    // The two lines after this contain atoms information: atoms and atom types
    fgets(line, buffer, readFilePointer);
    int numberOfAtoms, atomTypes;
    sscanf(line, "%i %*s\n", &numberOfAtoms);
    numberOfPolymers = numberOfAtoms / numberOfMonomers;
    fprintf(writeFilePointer, "%i atoms\n", numberOfAtoms);
    fgets(line, buffer, readFilePointer);
    sscanf(line, "%i %*s\n", &atomTypes);
    fprintf(writeFilePointer, "%i atom types\n", atomTypes);

    InitializeOffsetAngles(areAnglesFound, anglesString, " "); // Initializing the offset angles

    // handling the bond information:
    //reading number of bonds:
    int bonds;
    fgets(line, buffer, readFilePointer);
    sscanf(line, "%i %*s\n", &bonds);
    bonds = bonds + numberOfPolymers * archCrossLinks.numberOfCrossLinks; // updating the bonds: adding the crosslinks
    fprintf(writeFilePointer, "%i bonds\n", bonds);
    //reading bond types:
    fgets(line, buffer, readFilePointer); // skipping the line
    int numberOfBondTypes = 1;
    if(useDifferentBondTypeForCrossLinks)
        numberOfBondTypes = 3;
    fprintf(writeFilePointer, "%i bond types\n", numberOfBondTypes);

    //copying all lines until the atoms section is reached:
    char* query = "Atoms";
    CopyLinesTillQuery(readFilePointer, writeFilePointer, query, false);
    //Atoms section reached; title not printed
}

//Handle the copying of the Atoms section
void CopyAtomsSection(FILE* readFilePointer, FILE* writeFilePointer)
{
    // Printing the title
    fprintf(writeFilePointer, "Atoms\n");
    // copying all lines:
    char* query = "Bonds";
    CopyLinesTillQuery(readFilePointer, writeFilePointer, query, false);
    // Bonds section reached; title not printed
}

//Handle the copying of the Velocities section; Not present for the initial configuration files
void CopyVelocitiesSection(FILE* readFIlePointer, FILE* writeFilePointer, bool skipSection)
{
    if(!skipSection)
        fprintf(writeFilePointer, "Velocities\n"); // writing title of the section

    char* query = "Bonds";
    CopyLinesTillQuery(readFIlePointer, writeFilePointer, query, skipSection);
}

//Handle the copying of the Bonds section:
void CopyBondsSection(FILE* readFilePointer, FILE* writeFilePointer)
{
    // Printing the title:
    fprintf(writeFilePointer, "Bonds\n");

    // reading till the end of file:
    int buffer = 200;
    char line[buffer];
    int bondCounter = 0; // counter to keep track of serial number of the bonds;

    while(true) 
    {
        fgets(line, buffer, readFilePointer);

        if(feof(readFilePointer) != 0) // feof returns 1 or some non zero value when EOF is reached; it is updated only after reading a line is tried when no new lines exist
            break;

        //checking if line is empty:
        if(line[0] == '\n')
        {
            fprintf(writeFilePointer, "%s", line);
            continue;
        }

        // line containing bond information:
        int bondType;
        int atom1, atom2; // the two bonded atoms' IDs for a particular bond
        sscanf(line, "%*i %i %i %i\n", &bondType, &atom1, &atom2);
        //debugging:
        // if(bondCounter == 399 || bondCounter == 400)
        //     printf("Bond Type = %i\n", bondType);
        if(bondType == 1)
            fprintf(writeFilePointer, "%i %i %i %i\n", ++bondCounter, bondType, atom1, atom2); // writing the line only if the bond type is 1

    }
    // reading the file finished
    //adding the extra cross link bonds for both rings:
    int crossLinkBondType = 1;
    if(useDifferentBondTypeForCrossLinks)
        crossLinkBondType = 3;
    for(int i = 0; i < archCrossLinks.numberOfCrossLinks; i++)
    {
        for(int n = 0; n < numberOfPolymers; n++)
            fprintf(writeFilePointer, "%i %i %i %i\n", ++bondCounter, crossLinkBondType, archCrossLinks.crossLinks[i][0] + n * numberOfMonomers, archCrossLinks.crossLinks[i][1] + n * numberOfMonomers); // for the nth ring polymer
    }
    //done!
}


//Reads the initial state file from the readFilePath and writes the new modified state file to the writeFilePath
void ReadAndWriteInitialState(const char* readFilePath, const char* writeFilePath, const bool skipVelocitiesSection)
{

    // opening files:
    FILE* readFilePointer; 
    if((readFilePointer = fopen(readFilePath, "r")) == NULL)
    {
        printf("Read File could not be found at %s!\n", readFilePath);
        exit(1);
    }
    FILE* writeFilePointer;
    if((writeFilePointer = fopen(writeFilePath, "w")) == NULL)
    {
        printf("Write File could not be opened at %s!\n", writeFilePath);
        exit(1);
    }

    // handling the header part of the file: before any major sections:
    CopyHeader(readFilePointer, writeFilePointer);
    // Atoms section:
    CopyAtomsSection(readFilePointer, writeFilePointer);
    //Bonds Section:
    CopyBondsSection(readFilePointer, writeFilePointer);

    //closing file:
    fclose(readFilePointer);
    fclose(writeFilePointer);
}

//compares two strings and returns their difference: debugging
int CompareStrings(const char* string1, const char* string2)
{
    return strcmp(string1, string2);
}

//Reads the initial state file and generates the initial state for a particular architecture
void GenerateInitialState(void)
{
    char* readFilePath;
    char* writeFilePath;
    SetReadFilePath(&readFilePath);
    SetWriteFilePath(&writeFilePath);

    ReadAndWriteInitialState(readFilePath, writeFilePath, false);
    printf("Data file generated for %s at %s\n", architecture, writeFilePath);
    free(readFilePath);
    free(writeFilePath);
}

// Frees the allocated memoery for various global variables
void FreeMemory(void)
{
    free(offsetAngles);
    free(directory);
}

int main(int argc, char** argv)
{
    SetConstants(argc, argv);
    InitializeArchitecture(&archCrossLinks, architecture); // Initializes the Architecture struct to store the cross links information
    char* databaseFilePath;
    SetCrossLinksDatabaseFilePath(&databaseFilePath, numberOfMonomers); // The 
    ReadCrossLinks(&archCrossLinks, databaseFilePath); // Reads the cross links for the architecture from the database
    free(databaseFilePath);
    // printf("Number of Cross Links in %s : %i\n", architecture, archCrossLinks.numberOfCrossLinks);
    // printf("One example index: %i\n", archCrossLinks.crossLinks[archCrossLinks.numberOfCrossLinks - 1][1]);

    GenerateInitialState(); // Reads the Arc-0 initial state file and writes the new initial state file for the particular architecture
    FreeMemory();
    FreeArchitecture(&archCrossLinks);
}
