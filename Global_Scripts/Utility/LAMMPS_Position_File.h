#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#ifndef LAMMPS_POSITION_FILE_H // include guard
    #define LAMMPS_POSITION_FILE_H
    // This is a header file that contains objects and methods useful while reading the position data dumped during a simulation by the LAMMPS software.
    // The format of the dump file is the following:
    // - The dump file consists of "entries" at different timesteps of the simulation. Each entry has a few header lines that speicify certain details 
    //    of the simulation at that state followed by a whitespace separated table of the data of different atoms
    // - The header consists of:
        //ITEM: TIMESTEP
        // 10000
        // ITEM: NUMBER OF ATOMS
        // 200
        // ITEM: BOX BOUNDS ff ff ff
        // 0.0000000000000000e+00 4.7999999999999998e+00
        // 0.0000000000000000e+00 4.7999999999999998e+00
        // 0.0000000000000000e+00 2.0800000000000001e+01
    // This is followed by a header to the data table and then the actual data

    #define MAX_POLYMERS 2 // The maximum number of polymers expected in a simulation data file being read


    typedef struct atom_location atom_location;

    // variables necessary to specify a positions data file:
    struct simulation_data { // Common properties of the system at a particular timestep
        int timeStep; // The iteration number of the simulation at which the read has been read
        int numberOfMonomers; // The total number of monomers/atoms in the simulation; OR the maximum atomID present in the simulation
        int readNumberOfMonomers; // The total number of monomers present at a particular snapshot; read from the positions file
        double** bead_positions; // A look up table to store the bead positions of each atom belonging to different atom types
        // first index: atom ID-1
        // second index: x, y or z position
        int* atomTypes; // An array to store the type of each atom according to their atom ID; the index is atomID-1
        // aka the location in the bead_positions matrix for each atom in the simulation. This will be useful to track an atom by their ID.
        char* readFilePath; // A string to store the file path of the position file to be read
        FILE* readFilePointer; // The file pointer to the position file being read

        double box_bounds[3][2]; // The edge bounds of the simulation box in the x, y, and z directions
        int printingInterval; // The interval of iterations at which the position data is printed to the file
        bool isPrintingIntervalSet; // A flag to indicate whether the printing interval has been calculated on reading the position file

        bool suppressWarning; // A flag to indicate whether the warning for unequal number of monomers should be suppressed
    };
    // Note: accessing members of a struct:
    // struct simulation dataObject: dataObject.member
    // struct simulation dataPointer: dataPointer->member


    // type definition to use simulation_data as a variable type
    typedef struct simulation_data simulation_data;

    // Useful Methods to read a file:

    // Allocate memory to an array with the size passed
    static void AllocateArray(int** pointerToArray, int arrayLength)
    {
        *pointerToArray = (int*) calloc(arrayLength, sizeof(**pointerToArray));
    }

    // Allocates memory to the bead_positions look up table. Needs to be done only once per data file
    static void AllocateBeadPositions(simulation_data *data)
    {
        data->bead_positions = (double **) malloc(data->numberOfMonomers * sizeof(double*));
        assert(data->bead_positions); // checking if the memoery was allocated correctly. If not, throws an error
        for(int n = 0; n < data->numberOfMonomers; n++)
        {
            data->bead_positions[n] = (double*) malloc(3 * sizeof(double));
            assert(data->bead_positions[n]);
        }
        // alternative to doing this: use memcpy? Does that only work for arrays (and not look up tables)?
    }

    // Initialized the variables of the struct that are required to initialize the rest of the struct and subsequently, to read the position data properly
    static void InitializeRequiredSimulationDataFields(simulation_data *data, int numberOfMonomers, char* readFilePath)
    {
        // setting the numberOfMonomers
        data->numberOfMonomers = numberOfMonomers;
        
        // allocating memory for the readFilePath variable:
        size_t stringSize = strlen(readFilePath) + 1; // +1 to include the terminating null character as well
        data->readFilePath = (char*) malloc(stringSize * sizeof(char));
        // copying the readFilePath:
        strcpy(data->readFilePath, readFilePath);
    }

    // Initializes the derived quantities of the struct simulation_data; assuming the prerequisite quantities have already been set
    static void InitializeDerivedQuantities(simulation_data *data)
    {
        // timesteps information:
        data->timeStep = -1;
        data->isPrintingIntervalSet = false;
        data->printingInterval = 0;
        // read number of monomers:
        data->readNumberOfMonomers = data->numberOfMonomers; // unset
        // Warning message flag:
        data->suppressWarning = false;
        // initializing the atomTypes array:
        AllocateArray(&(data->atomTypes), data->numberOfMonomers);

        data->readFilePointer = fopen(data->readFilePath, "r");
        if(data->readFilePointer == NULL)
        {
            printf("ERROR: The position data file %s could not be read!\n", data->readFilePath);
            exit(1);
        }
        AllocateBeadPositions(data);
        // box bounds should already be initialized at declaration
    }

    // Initializes the simulation data struct completely
    void InitializeSimulationData(simulation_data *data, int numberOfMonomers, char* readFilePath)
    {
        InitializeRequiredSimulationDataFields(data, numberOfMonomers, readFilePath);
        InitializeDerivedQuantities(data);
    }

    // Frees the allocated memoery for the bead positions look up table
    static void FreeBeadPositions(simulation_data *data)
    {   
        for(int n = 0; n < data->numberOfMonomers; n++)
        {
            free(data->bead_positions[n]);
        }
        free(data->bead_positions);
    }

    // Frees the allocated memory for different members of the struct
    void FreeAllocatedMemory(simulation_data *data)
    {
        FreeBeadPositions(data);
        free(data->readFilePath);
        free(data->atomTypes);
    }

    //Skip the position data for a single timestep while reading the positions file
    void SkipPositions(simulation_data *data)
    {
        //assuming the reading cursor is at the start of a new timestep entry in the file:
        // skipping 9 lines of the header to get to the position data:
        int maxLineSize = 100;
        char buffer[maxLineSize];
        for(int i = 0; i < 9; i++)
        {
            fgets(buffer, maxLineSize, data->readFilePointer); // assuming each line has a return within maxLineSize number of characters
        }

        //skipping coord info:
        for(int i = 0; i < data->numberOfMonomers; i++)
        {
            fgets(buffer, maxLineSize, data->readFilePointer); // assuming each line has a return within maxLineSize number of characters
        }

    }

    // Initializes/Resets the values of an integer array
    static void ResetArray(size_t size, int* array)
    {
        for(int i = 0; i < size; i++)
            array[i] = 0;
    }


    // Reads the header and data for a single timestep and stores the data in the passed simulation_data object; returns true as long as EOF is not reached
    bool ReadPositions(simulation_data *data)
    {
        //assuming the reading cursor is at the start of a new timestep entry in the file:
        // skipping the first header line: "ITEM: TIMESTEP"
        int maxLineSize = 100;
        char buffer[maxLineSize];
        fgets(buffer, maxLineSize, data->readFilePointer);
        //Reading the timestep:
        fgets(buffer, maxLineSize, data->readFilePointer);
        int timeStep = 0;
        int successes = sscanf(buffer, "%i\n", &timeStep);

        // Setting the value of the printing interval if it is not already set:
        if(!data->isPrintingIntervalSet && data->timeStep != -1) // if the printing interval is unset and this isn't the first timestep being read
        {
            data->printingInterval = timeStep - data->timeStep; // Current timestep - previous timestep
            data->isPrintingIntervalSet = true;
        }
        
        // debugging:
        // printf("Timestep reached: %i\n", timeStep);
        
        // skipping the next header line: "ITEM: NUMBER OF ATOMS"
        fgets(buffer, maxLineSize, data->readFilePointer);
        // reading the number of atoms:
        fgets(buffer, maxLineSize, data->readFilePointer);
        int numberOfMonomers;
        successes = sscanf(buffer, "%i\n", &numberOfMonomers);

        // skipping the next header line: "ITEM: BOX BOUNDS"
        fgets(buffer, maxLineSize, data->readFilePointer);
        // reading the box bounds:
        double box_bounds[3][2];
        for(int i = 0; i < 3; i++)
        {
            fgets(buffer, maxLineSize, data->readFilePointer); // assuming each line has a return within maxLineSize number of characters
            successes = sscanf(buffer, "%lf %lf\n", &(box_bounds[i][0]), &(box_bounds[i][1]));
        }

        if(feof(data->readFilePointer))
        {
            printf("EOF reached\n");
            return false; // EOF reached
        }

        // Assigning the read values after checking that EOF is not reached: the values are not gibberish
        data->timeStep = timeStep;
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 2; j++)
                data->box_bounds[i][j] = box_bounds[i][j];
        }

        // Checking if the read number of monomers is as per expectation 
        if((numberOfMonomers != data->numberOfMonomers) && !data->suppressWarning)
        {
            data->suppressWarning = true;
            printf("WARNING: The read number of atoms %i is different from the one provided in the simulation data %i!\n", numberOfMonomers, data->numberOfMonomers);
            printf("This may be intentional in cases when the maximum atomID exceeds the number of atoms in the simulation.\n");
            printf("Continuing with the reading, but it should be noted that not all atoms in the bead_positions may be updated!\n");
            printf("Updating the read number of monomers in the simulation data based on the read value.\n");

            data->readNumberOfMonomers = numberOfMonomers;
        }
        
        // reading positions:
        int type;
        int atomID;
        double x, y, z; // variables to read the positions of a monomer

        // skipping the data header: ITEM: ATOMS id type x y z 
        fgets(buffer, maxLineSize, data->readFilePointer);
        
        for(int i = 0; i < numberOfMonomers; i++)
        {
            fgets(buffer, maxLineSize, data->readFilePointer);
            sscanf(buffer, "%i %i %lf %lf %lf\n", &atomID, &type, &x, &y, &z);
            // Checking validity of monomer atom ID:
            if(atomID > data->numberOfMonomers)
            {
                printf("ERROR: The quantity of atoms has exceeded the set value %i!\n", data->numberOfMonomers);
                exit(1);
            }

            // Storing the positions of the atoms:
            int atomIndex = atomID - 1;
            data->bead_positions[atomIndex][0] = x;
            data->bead_positions[atomIndex][1] = y;
            data->bead_positions[atomIndex][2] = z;

            // Storing the atom type:
            data->atomTypes[atomIndex] = type;
        }

        //the cursor should ideally be at the start of the next timestep data or at end of the file
        // EOF not reached
        return true;
    }

    // Closes the read file:
    void CloseReadFile(simulation_data *data)
    {
        fclose(data->readFilePointer);
    }

    // Writes the header and the data stored in the passed instance of the data struct (for a single timestep) to a file pointed to by a file pointer
    void WritePositions(simulation_data* data, FILE* writeFilePointer)
    {
        // Writing the header to the file:
        fprintf(writeFilePointer, "ITEM: TIMESTEP\n%i\n", data->timeStep); // Writing the timestep
        fprintf(writeFilePointer, "ITEM: NUMBER OF ATOMS\n%i\n", data->numberOfMonomers); // Writing the number of monomers/atoms
        fprintf(writeFilePointer, "ITEM: BOX BOUNDS ff ff ff\n"); // Currently only supports the fixed boundary condition
        for(int i = 0; i < 3; i++)
        {
            fprintf(writeFilePointer, "%lf %lf\n", data->box_bounds[i][0], data->box_bounds[i][1]); // Writing the box bounds
        }
        fprintf(writeFilePointer, "ITEM: ATOMS id type x y z\n"); // Writing the atoms positions header line
        // Writing the positions data:
        for(int i = 0; i < data->numberOfMonomers; i++)
        {
            fprintf(writeFilePointer, "%i %i %lf %lf %lf\n", i+1, data->atomTypes[i], data->bead_positions[i][0], data->bead_positions[i][1], data->bead_positions[i][2]);
        }
    }

    // Calculates the number of polymers present in the system by reading the first snapshot
    int ReadNumberOfPolymers(int numberOfMonomers, char* readFilePath)
    {
        simulation_data read_data;
        InitializeSimulationData(&read_data, MAX_POLYMERS * numberOfMonomers, readFilePath);

        // Reading first snapshot:
        ReadPositions(&read_data);
        int numberOfPolymers = read_data.readNumberOfMonomers / numberOfMonomers; // setting number of polymers

        // Freeing up memory for read_data struct:
        FreeAllocatedMemory(&read_data);
        return numberOfPolymers;
    }
#endif