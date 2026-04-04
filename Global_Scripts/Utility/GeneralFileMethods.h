#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

// including the system_file_paths.h file:
#include "../System_File_Paths/system_file_paths.h"

#ifndef GENERAL_FILE_METHODS_H // include guard
    #define GENERAL_FILE_METHODS_H
    //Append three strings together in order and puts the result in the pointed string location
    void Append(char** finalStringPointer, const char* string1, const char* string2, const char* string3)
    {
        int bytes = asprintf(finalStringPointer, "%s%s%s", string1, string2, string3);
        if(bytes == -1)
        {
            printf("asprintf failed to allocate memory to append strings.\n");
            exit(1);
        }
    }

    //Sets the folder name according to the run number. For example: "noConfinement/run11/"
    void SetFolderName(char** folderNamePointer, char* folderNamePrefix, int runNumber)
    {   
        int bytes = asprintf(folderNamePointer, "%s/run%i/", folderNamePrefix, runNumber);
        if(bytes == -1)
        {
            printf("Memory could not be allocated for the folder name.\n");
            exit(1);
        }
    }

    //Sets the folder name according to the run number and uncut flag. For example: "noConfinement/run11_uncut/"
    void SetUncutFolderName(char** folderNamePointer, char* folderNamePrefix, int runNumber, bool isSystemUncut)
    {
        char* uncutLabel = "";
        if(isSystemUncut)
            uncutLabel = "_uncut";
        int bytes = asprintf(folderNamePointer, "%s/run%i%s/", folderNamePrefix, runNumber, uncutLabel);
        if(bytes == -1)
        {
            printf("Memory could not be allocated for the folder name.\n");
            exit(1);
        }
    }

    //Sets the directory path by taking into account the number of Monomers and the special simulation variable set in the sysPaths file.
    // For example Adds "b500/" at the end of the passed path or goes to a simulation folder under b500/Previous_Attempts/
    void SetDirectoryPath(char** directoryPathPointer, char* pathPrefix, int numberOfMonomers, char* specialSimulation)
    {
        int bytes;
        if(strlen(specialSimulation) == 0) // SPECIAL_SIMULATION not set; default location
        {
            bytes = asprintf(directoryPathPointer, "%sb%i/", pathPrefix, numberOfMonomers);
        }
        else // SPECIAL_SIMULATION set
        {
            bytes = asprintf(directoryPathPointer, "%sb%i/%s/", pathPrefix, numberOfMonomers, specialSimulation);
        }
        if(bytes == -1)
        {
            printf("Memory could not be allocated for the directory path!\n");
            exit(1);
        }
    }

    // Sets the directory path in the Rg folder by taking into account the number of Monomers and the Rg simulation variable set in the sysPaths file.
    // For example Adds "b500/" at the end of the passed path or goes to a simulation folder under b500/Previous_Attempts/
    void SetRgDirectoryPath(char** directoryPathPointer, char* pathPrefix, int numberOfMonomers)
    {
        int bytes;
        if(strlen(R_G_SIMULATION) == 0) // SPECIAL_SIMULATION not set; default location
        {
            bytes = asprintf(directoryPathPointer, "%sb%i/", pathPrefix, numberOfMonomers);
        }
        else // SPECIAL_SIMULATION set
        {
            bytes = asprintf(directoryPathPointer, "%sb%i/%s/", pathPrefix, numberOfMonomers, R_G_SIMULATION);
        }
        if(bytes == -1)
        {
            printf("Memory could not be allocated for the directory path!\n");
            exit(1);
        }
    }

    // Sets the read file path for the position dump of the simulation
    // Args:
    // - filePathPointer: The pointer to the string where the file path will be stored (does not need to be allocated)
    // - directory: The simulation directory path where the simulation files are, up to the b<N>/ directory; where N is the number of monomers
    // - architecture: The name of the architecture used in the simulations
    // - runIndex: The index to choose which of the independent runs to read from
    // - forInitialization: A flag to indicate whether the simulations run were for initializing the mixed state of polymers
    void SetPositionDumpFilePath(char** filePathPointer, char* directory, char* architecture, int runIndex, bool forInitialization)
    {
        char* fileName;
        if(forInitialization)
            fileName = "distribution_positions.dump";
        else
            fileName = "visual.dump";

        char* folder; 
        SetFolderName(&folder, architecture, runIndex);
        // Setting the final path:
        Append(filePathPointer, directory, folder, fileName);
        free(folder);
        if(*filePathPointer == NULL)
        {
            printf("ERROR: Memory could not be allocated for the position dump file path.\n");
            exit(1);
        }
    }

    // Sets the file path for a custom file to write data to
    // Args:
    // - filePathPointer: The pointer to the string where the file path will be stored (does not need to be allocated)
    // - directory: The simulation directory path where the simulation files are, up to the b<N>/ directory; where N is the number of monomers
    // - architecture: The name of the architecture used in the simulations
    // - runIndex: The index to choose which of the independent runs to read from
    // - fileName: The name of the file to be written to with the extension
    void SetCustomWriteFilePath(char** filePathPointer, char* directory, char* architecture, int runIndex, char* fileName)
    {
        char* folder; 
        SetFolderName(&folder, architecture, runIndex);
        // Setting the final path:
        Append(filePathPointer, directory, folder, fileName);
        free(folder);
        if(*filePathPointer == NULL)
        {
            printf("ERROR: Memory could not be allocated for the custom write file path.\n");
            exit(1);
        }
    }

    //Reads the optional arguments (destination architecture and runIndex) that are passed after the mandatory ones and stores them in the first parameter passed to this function. Returns the runIndex 
    int GetOptionalArguments(char** destinationArchitecture, int argc, char** argv, const int numberOfMandatoryArguments)
    {
        if(argc <= numberOfMandatoryArguments + 1) // checking if only the mandatory arguments are passed
        {
            printf("No optional arguments passed.\n");
            return -1;
        }
        else if(argc == numberOfMandatoryArguments + 1 + 1) // 1 mandatory argument passed
        {
            char* argument = argv[numberOfMandatoryArguments + 1];
            //Checking if the argument passed is an integer, and thus, the runIndex:
            int runIndex = atoi(argument);
            if(runIndex == 0)
            {
                //argument is not a number
                printf("No argument for the run index passed. Continuing without specifying a particular run\n");
                strcpy(*destinationArchitecture, argument);
                return -1;
            }
            else
            {
                printf("No argument for the destination architecture passed. Retaining the same architecture as read from mixed state file\n");
                return runIndex;
            }
        }
        else // both arguments passed
        {
            //assuming the aarchitecture is passed first and then the runIndex
            strcpy(*destinationArchitecture, argv[numberOfMandatoryArguments + 1]);
            char* runIndexArg = argv[numberOfMandatoryArguments + 1 + 1];

            int runIndex = atoi(runIndexArg);
            if(runIndex == 0)
            {
                printf("The argument for the run Index \"%s\" could not be converted to an integer. Please enter a valid number\n", runIndexArg);
                return -1;
            }
            else
                return runIndex;

        }
    }

    //Reads the optional runIndex in the arguments array passed from the command line and returns it
    int GetRunIndex(int argc, char** argv, const int numberOfMandatoryArguments)
    {
        if(argc <= numberOfMandatoryArguments + 1) // checking if only the mandatory arguments are passed
        {
            printf("No argument for the run index passed. Continuing without specifying a particular run\n");
            return -1;
        }
        else // run Index passed
        {
            int runIndex = atoi(argv[numberOfMandatoryArguments + 1]); // returns 0 if the argument cannot be converted to an integer
            
            if(runIndex == 0)
            {
                printf("The argument for the run Index \"%s\" could not be converted to an integer. Please enter a valid number\n", argv[numberOfMandatoryArguments + 1]);
                exit(1);
            }
            else
                return runIndex;
        }
    }

    //Sets the position file path based on the passed parameters: Unused
    void SetPositionDumpPath(char** filePathPointer ,int numberOfMonomers, char* architecture, int runNumber, char* positionDumpDirectory)
    {
        char* fileName;
        if(numberOfMonomers == 200)
            fileName = "positions.dump";
        else if(numberOfMonomers == 500)
            fileName = "visual.dump";

        char* folderName;

        SetFolderName(&folderName, architecture, runNumber);
        Append(filePathPointer, positionDumpDirectory, folderName, fileName);
        free(folderName);
    }

    // Sets the file path to the cross-links database file based on the number of monomers passed
    void SetCrossLinksDatabaseFilePath(char** filePathPointer, int numberOfMonomers)
    {
        int bytes = asprintf(filePathPointer, "%sb%i/Arc_Cross_Links_%i.txt", HOME_DIR CREATE_INITIAL_STATES, numberOfMonomers, numberOfMonomers);
        if(bytes == -1)
        {
            printf("Memory could not be allocated for the cross-links database file name!\n");
            exit(1);
        }
    }

    // Sets the file path to the diamaters database file based on the number of monomers passed
    void SetDiametersDatabaseFilePath(char** filePathPointer, int numberOfMonomers)
    {
        int bytes = asprintf(filePathPointer, "%sb%i/Diameters.csv", HOME_DIR SEGREGATION, numberOfMonomers);
        if(bytes == -1)
        {
            printf("Memory could not be allocated for the diameters database file name!\n");
            exit(1);
        }
    }

    // Sets the file path to the axis length database file based on the number of monomers passed
    void SetAxisLengthDatabaseFilePath(char** filePathPointer, int numberOfMonomers)
    {
        int bytes = asprintf(filePathPointer, "%sb%i/AxisLengths.csv", HOME_DIR SEGREGATION, numberOfMonomers);
        if(bytes == -1)
        {
            printf("Memory could not be allocated for the axis length database file name!\n");
            exit(1);
        }
    }

    // Returns whether the SPECIAL_SIMULATION specified in system_file_paths.h involves an infinite cylinder for confinement
    bool IsCylinderInfinite(void)
    {
        char query[] = "inf_";
        return strncmp(SPECIAL_SIMULATION, query, strlen(query)) == 0; // If 'inf_' exists as the start of the special simulation string
    }

    // Prints a string array
    void PrintArray(size_t arraySize, char** stringArray)
    {
        for(int i = 0; i < arraySize - 1; i++)
        {
            printf("%s, ", stringArray[i]);
        }
        // Printing last string separately:
        printf("%s\n", stringArray[arraySize - 1]);
    }

    // Prints the argument instructions for each directory keyword
    void PrintKeywordInstructions(int numberOfKeywords, char** acceptedKeywords, char** keywordInstructions)
    {
        printf("Instructions for the additional arguments for each of the keywords:\n");
        for(int i = 0; i < numberOfKeywords; i++)
        {
            printf("%s: ", acceptedKeywords[i]);
            printf("%s\n", keywordInstructions[i]);
        }
    }

    // Searches for a query string in an array of string and returns the fouund index. If the query is not found, then returns -1
    // Checks for equality of strings by ignoring case
    int SearchStringArray(char* query, size_t arraySize, char** stringArray)
    {
        int foundIndex = -1;
        for(int i = 0; i < arraySize; i++)
        {
            if(strcasecmp(query, stringArray[i]) == 0) // Match found
            {
                return i;
            }
        }
        return foundIndex;
    }

    // Prints the formatted string to the stdout and to the passed file stream
    void bothprintf(FILE* writeFilePointer, char const *fmt, ...)
    {
        va_list ap; // defining the pointer that will travel the variable argument list fmt

        // Printing to stdout
        va_start(ap, fmt); // initialises the pointer to the start of the variable argument list fmt
        vprintf(fmt, ap);
        va_end(ap); // Cleans up the pointer so that ap cannot be used unless initialised
        // Printing to file:
        va_start(ap, fmt); // initialises the pointer to the start of the variable argument list fmt
        vfprintf(writeFilePointer, fmt, ap);
        va_end(ap);
    }
#endif