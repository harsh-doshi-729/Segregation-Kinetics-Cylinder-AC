#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#ifndef SEGREGATION_TIMES_FILE_H // include guard
    #define SEGREGATION_TIMES_FILE_H
    #include "../System_File_Paths/system_file_paths.h"
    #include GENERAL_FILE_METHODS

    #define CRITERIA_LENGTH 14

    // This header file contains objects and mehtods useful for reading the segregation times for multiple runs of a single architecture.
    // The segregation times file is expected to have the following format:
    // Run Index, Segregation Time steps (f##_s##_t###)
    // 1, #####
    // 2, #####
    // and so on...

    struct segregation_times {
        char* filePath; // The file path to the segregation times file
        bool isFilePathInternal; // A flag to indicate whether the file path has been set my methods of this header file or externally
        int* segregationTimes; // An array to store the segregation times indexed with the run number
        int numberOfRuns; // The number of runs stored in the file

        char* segregationCriteria; // The criteria used for segregation of the two polymers
    };

    typedef struct segregation_times segregation_times;

    // Allocates the memory for the segregation times array
    void AllocateSegregationArray(segregation_times* data)
    {
        data->segregationTimes = (int*) calloc(data->numberOfRuns, sizeof(*(data->segregationTimes)));
        // setting default value of segregation time to -1 (did not segregate)
        for(int i = 0; i < data->numberOfRuns; i++)
        {
            data->segregationTimes[i] = -1; 
        }
    }

    // Sets the filepath to the segregation times csv file
    void SetSegregationTimesFilePath(char** filePathPointer, int numberOfMonomers, char* architecture, char* segregationCriterion)
    {
        char* directory;
        SetDirectoryPath(&directory, NEW_SEGREGATION, numberOfMonomers);
        int bytes = asprintf(filePathPointer, "%s%s/Analysis/segregationTimes_%s.csv", directory, architecture, segregationCriterion);
        if(bytes == -1)
        {
            printf("ERROR: Memory could not be allocated for the segregation times file path! Terminating.\n");
            exit(1);
        }
    }


    // Reads the segregation times file and stores the necessary information in the struct
    // Reads the segregation times at the passed filePath and sets them in the passed integer array
    void ReadSegregationTimes(segregation_times* data)
    {
        FILE* segregationFile = fopen(data->filePath, "r");
        if(segregationFile == NULL)
        {
            printf("ERROR: The file %s could not be opened! Terminating.\n", data->filePath);
            exit(1);
        }

        // Reading file:
        int buffer = 200; // Length of the string buffer for lines to be read from the file
        char line[buffer]; // The buffer into which a line will be read

        // reading segregation criteria from the header:
        fgets(line, buffer, segregationFile);
        int successes = sscanf(line, "Run Index, Segregation Time steps (%[^()])\n", data->segregationCriteria);
        // The format specifier %[^chars] reads all characters except for the ones specified after the ^
        if(successes != 1)
            printf("WARNING: The segregation criteria could not be read from the segregation times file!\n");

        // Reading segregation times from each line
        while(true)
        {
            // Reading line:
            fgets(line, buffer, segregationFile);

            // exit condition: EOF reached
            if(feof(segregationFile))
            {
                printf("Finished reading Segregation Times file!\n");
                break;
            }

            // Reading segregation time from buffer line:
            int segregationTime;
            int runNumber;
            int successes = sscanf(line, "%i, %i\n", &runNumber, &segregationTime);
            if(successes != 2)
            {
                printf("WARNING: An unexpected number of arguments were read from a line in segregationTimes.csv.\n");
            }
            if(runNumber > data->numberOfRuns)
            {
                printf("ERROR: Run %i encountered in the segregation times file exceeded the set number of runs of %i! Terminating.\n", runNumber, data->numberOfRuns);
                exit(1);
            }
            data->segregationTimes[runNumber - 1] = segregationTime; // setting the value in the array
            
        }

    }

    // Initializes and reads the segregation_times struct by assuming the default location for the segregation file:
    void InitializeSegregationTimes(segregation_times* data, int numberOfMonomers, char* architecture, int numberOfRuns)
    {
        data->isFilePathInternal = true;
        data->numberOfRuns = numberOfRuns;
        
        // Allocating memory for the segregation criteria string
        data->segregationCriteria = (char*) malloc(CRITERIA_LENGTH * sizeof(char));
        // Setting Segregation criteria using the defined macro in system_file_paths.h
        strcpy(data->segregationCriteria, SEGREGATION_CRITERION);
        SetSegregationTimesFilePath(&(data->filePath), numberOfMonomers, architecture, data->segregationCriteria);

        AllocateSegregationArray(data);
        ReadSegregationTimes(data);
    }

    // Initialiizes the segregation_times struct by accepting the absolute path to the segregation file
    void InitializeSegregationTImesWithFilePath(segregation_times* data, char* filePath, int numberOfRuns)
    {
        data->filePath = filePath;
        data->isFilePathInternal = false;
        data->numberOfRuns = numberOfRuns;

        AllocateSegregationArray(data);
        ReadSegregationTimes(data);
    }

    // Prints the segregation times to console: for debugging
    void PrintSegregationTimes(segregation_times* data)
    {
        printf("Segregation Times:\n");
        for(int i = 0; i < data->numberOfRuns; i++)
        {
            printf("Run %i:\t%i\n", i+1, data->segregationTimes[i]);
        }
    }

    // Frees the memory allocated for the member variables of the struct
    void FreeSegregationTimes(segregation_times* data)
    {
        free(data->segregationTimes);
        free(data->segregationCriteria);
        if(data->isFilePathInternal)
            free(data->filePath);
    }
#endif