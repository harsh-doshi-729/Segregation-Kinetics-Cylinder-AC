// the output file generated from the fix print command seems to be buggy. It prints for every time an MC step is attempted. And it also prints a few extra steps in the start

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "../../../../Scripts/System_File_Paths/system_file_paths.h"
#include GENERAL_FILE_METHODS

int Index;
char* architecture;
char* directory; //the directory where the architecture folders that contain the simulation files are located
char* readFilePath; //relative path
char* writeFilePath; //relative path
int offset; // skip the first few entries
int skipInterval; // skip these many entries after each entry is read and written
int numberOfMCSteps;
int run;
int numberOfMonomers;


//Sets the global constant variables based on the arguments passed
void SetConstants(int argc, char** argv)
{
    // parsing the arguments: expecting just 1 argument for the architecture
    if(argc < 3) // only the default program was passed
    {
        printf("Not enough arguments passed! Please pass the architecture(string) and number of monomers(int) as an argument from the command line\nFor example:\n./correctExecuatable.out Arc2 500\n");
        exit(1);
    }
    else
    {
        //assuming atleast two arguments were passed; aprat from the default program name
        architecture = argv[1]; // the first argument after the program name
        numberOfMonomers = atoi(argv[2]);
        if(numberOfMonomers == 0)
        {
            printf("The second argument %s could not be converted to an integer! Please pass a valid number for the number of monomers\n", argv[2]);
	        exit(1);
        }
    }

    Index = 1;
    SetDirectoryPath(&directory, NEW_SEGREGATION, numberOfMonomers);
    offset = 2; // skip the first few entries
    skipInterval = 2 * numberOfMonomers + 1; // skip these many entries after each entry is read and written; this should be equal to 1 less than number of MC moves attempted in one timestep
    numberOfMCSteps = 1000;
    run = 8;
}

//method to set the FilePath to be read/written
void SetFilePath(char** filePathPointer, const int run, const int index, const char* prefix)
{
    int bytes = asprintf(filePathPointer, "%s%s/run%i/%scom%i.dat", directory, architecture, run, prefix, index);
    if(bytes == -1)
    {
        printf("Memory allocation for asprintf failed\n");
        exit(1);
    }
}

//method to free the filePaths
void FreeFilePaths(void)
{
    free(readFilePath);
    free(writeFilePath);
}

// method to read from the file and write the corrected data
void CorrectOutput(void)
{
    FILE* readFilePointer = fopen(readFilePath, "r");
    FILE* writeFilePointer = fopen(writeFilePath, "w");

    int bufferSize = 200;
    char line[bufferSize];

    // writing the first header line as is:
    fgets(line, bufferSize, readFilePointer);
    fprintf(writeFilePointer, "%s", line);

    //skipping the next two lines
    for(int i = 0; i < offset; i++)
        fgets(line, bufferSize, readFilePointer);
    
    bool isMCRegion = true; // the file starts with MC

    int time;
    double x, y, z;
    //writing each line and skipping the next  lines:
    while(isMCRegion)
    {
        fscanf(readFilePointer, "%i %lf %lf %lf\n", &time, &x, &y, &z);
        // the line has been read and stored into line: printing to file
        //fprintf(writeFilePointer, "%s", line);
        fprintf(writeFilePointer, "%i %.12lf %.12lf %.12lf\n", time, x, y, z);

        //checking if the last step of the MC simulation is reached
        if(time == numberOfMCSteps)
        {
            isMCRegion = false;
            skipInterval = skipInterval - 2; // changing the skip interval so that the entries in the Langevin region are not skipped
        }
        //skipping lines:
        for(int i = 0; i < skipInterval; i++)
            fgets(line, bufferSize, readFilePointer);
    }
    skipInterval = skipInterval + 2; // restoring to original skip interval

    //skipping 1 intermediate header line
    fgets(line, bufferSize, readFilePointer);
    //copying the rest of the file as is:
    while(fgets(line, bufferSize, readFilePointer) != NULL)
    {
        // the line has been read and stored into line: printing to file
        fprintf(writeFilePointer, "%s", line);
    }
    //EOF reached
    fclose(readFilePointer);
    fclose(writeFilePointer);
}

int main(int argc, char** argv)
{
    SetConstants(argc, argv);
    for(run = 1; run <= 50; run++)
    {
	for(Index = 1; Index <= 2; Index++)
    	{
            //setting the read and write file paths:
            SetFilePath(&readFilePath, run, Index, "");
            SetFilePath(&writeFilePath, run, Index, "corrected_");
            //reading and writing the entries:
            CorrectOutput();
            //freeing filePaths:
            FreeFilePaths();
    	}
    printf("CoM Correction: Run %i done!\n", run);
    }
    return 0;
}
