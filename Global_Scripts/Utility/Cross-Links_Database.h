// A header file that helps with the reading of the cross-links data base file that contains the cross-link information for various architectures
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include MODULO

#ifndef CROSS_LINKS_DATABASE_H // include guard
    #define CROSS_LINKS_DATABASE_H
    #include "GeneralFileMethods.h" // present in the same directory

    #define LIN_ARC_STRING "Arc_Lin" // The string to search for in linear architectures
    #define ASPECT_RATIO 5 // The ratio of the axis length to the diameter of the confinment cylinder when the cylinder is finite

    struct Architecture {
        char* architectureName;
        int numberOfCrossLinks;
        // int numberOfMonomers; // The number of monomers present in a single polymer; TODO: Add this?
        int** crossLinks; // A variable length array of two-integer arrays to store the monomer indices being cross linked
        double confinementDiameter; // The diameter of the confining cylinder so that the R_g / D ratio stays roughly constant for all architectures
        double axisLength; // The axis length of the confining cylinder; 
        // in case of infinite cylinders, this is the effective span of the polymer with a particular architecture along the long axis
        bool isArcLinear; // A boolean variable that indicates whether the architecture is a linear polymer or based on a linear polymer
    };

    typedef struct Architecture Architecture;

    // Determines whether the architecture is a linear architecture based on the architecture name
    bool IsArchitectureLinear(char* architectureName)
    {
        if(strstr(architectureName, LIN_ARC_STRING) != NULL) // strstr returns a pointer to the first occurence of the string if it is found
            return true;
        else
            return false;
    }

    // Initializes the Architecture struct with the name of the architecture
    void InitializeArchitecture(Architecture* architecture, char* architectureName)
    {
        architecture->architectureName = architectureName;
        architecture->confinementDiameter = -1; // default value
        architecture->isArcLinear = IsArchitectureLinear(architectureName);
    }

    // Allocates memory for the crossLinks array
    void AllocateCrossLinks(Architecture *architecture)
    {
        architecture->crossLinks = (typeof(architecture->crossLinks)) malloc(architecture->numberOfCrossLinks * sizeof(*(architecture->crossLinks)));
        for(int i = 0; i < architecture->numberOfCrossLinks; i++)
            architecture->crossLinks[i] = (typeof(*(architecture->crossLinks))) malloc(2 * sizeof(**(architecture->crossLinks)));
    }

    // Reads the cross-links database file and fills the cross links in the struct
    void ReadCrossLinks(Architecture* architecture, char* dataFilePath)
    {
        // reading the file:
        FILE* filePointer;
        if((filePointer = fopen(dataFilePath, "r")) == NULL)
        {
            printf("ERROR: Cross Links Database file could not be opened!\nPath: %s\n", dataFilePath);
            exit(1);
        }

        bool hasArchBeenFound = false;

        int buffer = 500;
        char line[buffer];
        while(!hasArchBeenFound) // reading all lines till the architecture is found
        {
            fgets(line, buffer, filePointer);

            //checking if EOF is reached:
            if(feof(filePointer) != 0)
                break; //exiting the loop

            if(line[0] == '#' || line[0] == '\n') // comment line or blank line
                continue;
            else // must have encountered an architecture line
            {
                int noOfCrossLinks; // local copies of the architecture and crosslink variables

                char arch[buffer]; // reading the architecture in each architecture line
                sscanf(line, "%s %i", arch, &noOfCrossLinks); // parsing the line contents
                // printf("Read Architecture: %s\n", arch); // Printing the read architecture; for debugging

                if(strcmp(arch, architecture->architectureName) == 0) // architecture found
                {
                    // reading the cross links into the array:
                    architecture->numberOfCrossLinks = noOfCrossLinks;
                    AllocateCrossLinks(architecture); // allocating memory for cross links
                    for(int i = 0; i < noOfCrossLinks; i++)
                    {
                        fgets(line, buffer, filePointer);

                        if(feof(filePointer)) // Checking if the end of file is reached prematurely
                        {
                            printf("ERROR: End of the cross-links database file reached while parsing crosslink for architecture %s!", architecture->architectureName);
                            printf("Unexpected EOF encountered; please check the database file. Terminating.\n");
                            exit(1);
                        }

                        sscanf(line, "%i %i\n", &(architecture->crossLinks[i][0]), &(architecture->crossLinks[i][1])); // scanning monomer indices
                    }
                    hasArchBeenFound = true;
                }
                else // incorrect architecture:
                {
                    // skip the cross link lines:
                    for(int i = 0; i < noOfCrossLinks; i++)
                        fgets(line, buffer, filePointer);
                }

            }

        }

        fclose(filePointer);
        
        if(!hasArchBeenFound)
        {
            printf("Architecture %s was not found in the database!\n", architecture->architectureName);
            exit(1);
        }
        else
        {
            printf("Architecture %s was found in the database with %i cross links\n", architecture->architectureName, architecture->numberOfCrossLinks);
        }
    }

    // Returns the axis length of the cylinder in case the cylinder is 'infinite'
    double GetInfiniteCylinderAxisLength(int numberOfMonomers)
    {
        if(numberOfMonomers == 200)
            return 500; // Length of long cylinder in case of 200 monomers
        else if(numberOfMonomers == 500)
            return 1000; // Length of long cylinder in case of 500 monomers
    }

    // Reads the cylinder diamaters database file and finds the diameter corresponding to current architecture
    void ReadDiameter(Architecture* architecture, char* dataFilePath)
    {
        FILE* filePointer = fopen(dataFilePath, "r");
        if(filePointer == NULL)
        {
            printf("ERROR: The diamater database file could not be opened!\nPath: %s\n", dataFilePath);
            exit(1);
        }

        int buffer = 200; // Size of buffer to which a line is read
        char line[buffer];
        // skipping header:
        fgets(line, buffer, filePointer);

        // reading each architecture line:
        bool hasArchBeenFound = false;
        // Variables to store the read data:
        char arch[buffer]; 
        double diameter;

        while(!hasArchBeenFound)
        {
            fgets(line, buffer, filePointer);
            // Checking for EOF
            if(feof(filePointer) != 0) // EOF reached
            {
                break;
            }
            // Reading data using strtok:
            char* delimiters = ",\n"; // comma and a return character
            char* token = strtok(line, delimiters);
            strcpy(arch, token); // saving first token as architecture
            token = strtok(NULL, delimiters); // Skippng the second column
            token = strtok(NULL, delimiters);
            diameter = atof(token);
            
            if(strcmp(arch, architecture->architectureName) == 0) // Architecture found
            {
                hasArchBeenFound = true;
                architecture->confinementDiameter = diameter;
            }
        }
        if(hasArchBeenFound)
        {
            printf("Architecture %s was found in the diameters database with a diameter of %lf.\n", architecture->architectureName, architecture->confinementDiameter);
        }
        else
        {
            printf("ERROR: The architecture %s was not found in the diameters database.\n", architecture->architectureName);
            exit(1);
        }
    }

    // Reads the axis length (or effective axis length) of the confining cylinder corresponding to the architecture passed from the database
    // The value of the axis length is set in the Architecture struct instance.
    // Args:
    // - architecture: The pointer to the Architecture struct instance for which the data is to be read; assuming diameter has been set
    // - dataFilePath: The file path to the axis lengths database file listing the axis lengths for the combination of architetcure and confining diameter
    void ReadAxisLength(Architecture* architecture, char* dataFilePath)
    {
        // Opening database file to read:
        FILE* filePointer = fopen(dataFilePath, "r");
        if(filePointer == NULL)
        {
            printf("ERROR: The axis length database file could not be opened!\nPath: %s\n", dataFilePath);
            exit(1);
        }

        int buffer = 200; // Size of buffer to which a line is read
        char line[buffer];
        // skipping header:
        fgets(line, buffer, filePointer);

        // reading each architecture line:
        bool hasCombinationBeenFound = false; // A flag indicating whether the correct combination of architecture and diameter has been found in a single line
        // Variables to store the read data:
        char arch[buffer]; 
        double readDiameter; // The diameter that has been read from the axis length entry
        double axisLength;

        while(!hasCombinationBeenFound)
        {
            fgets(line, buffer, filePointer);
            // Checking for EOF
            if(feof(filePointer) != 0) // EOF reached
            {
                break;
            }
            // Reading data using strtok:
            char* delimiters = ",\n"; // comma and a return character
            char* token = strtok(line, delimiters);
            strcpy(arch, token); // saving first token as architecture
            token = strtok(NULL, delimiters); // Reading the second column: diameter
            readDiameter = atof(token);
            token = strtok(NULL, delimiters); // Reading the third column: axis length
            axisLength = atof(token);

            if(readDiameter == 0 || axisLength == 0)
            {
                printf("The diameter or axis lengths could not be read properly for '%s' from the axis length database file '%s'.\n", arch, dataFilePath);
                exit(1);
            }
            
            
            if((strcmp(arch, architecture->architectureName) == 0) && (architecture->confinementDiameter == readDiameter)) // Checking for exact equality of diameters since they should be exactly thr same
            {
                // Architecture and diameter combination found
                hasCombinationBeenFound = true;
                architecture->axisLength = axisLength;
            }
        }
        if(hasCombinationBeenFound)
        {
            printf("Architecture %s with diameter %lf was found in the axis length database with an axis length of %lf.\n", architecture->architectureName, architecture->confinementDiameter, architecture->axisLength);
        }
        else
        {
            printf("ERROR: The architecture %s and diameter %lf combination was not found in the axis length database.\n", architecture->architectureName, architecture->confinementDiameter);
            exit(1);
        }
    }


    // Reads the diameter and axis length of the confining cylinder corresponding to the architecture passed.
    // In case the simulation involves confinement in an infinite cylinder, the effective axis length is read from the database file
    // The effective axis length is the average length along the long axis of a cylinder spanned by a single polymer of the given architecture
    // If the cylinder is finite, the axis length is set based on the diameter and a predetermined aspect ratio
    // Sets the values of the diameter and axis length in the Architecture struct instance.
    // Args:
    // - architecture: The pointer to the Architecture struct instance for which the data is to be read
    // - numberOfMonomers: The number of monomers present in a single polymer of the system
    void ReadDiameterAndAxisLength(Architecture* architecture, int numberOfMonomers)
    {
        char* diameterDBFilePath;
        SetDiametersDatabaseFilePath(&diameterDBFilePath, numberOfMonomers);
        char* lengthDBFilePath;
        SetAxisLengthDatabaseFilePath(&lengthDBFilePath, numberOfMonomers);
        ReadDiameter(architecture, diameterDBFilePath);
        if(IsCylinderInfinite())
            ReadAxisLength(architecture, lengthDBFilePath);
        else
            architecture->axisLength = architecture->confinementDiameter * ASPECT_RATIO;

        free(diameterDBFilePath);
        free(lengthDBFilePath);
    }

    // Calculates the volume of the confining cylinder
    // In case of an infinite cylinder, the effective volume occupied by a single polymer is returned.
    // Args:
    // - architecture: The initialized Architecture struct (passed by value)
    double CalculateConfinementVolume(Architecture architecture)
    {
        // Volume = pi/4 * d^2 * l; d = diameter, l = axis length
        return M_PI / 4 * pow(architecture.confinementDiameter, 2) * architecture.axisLength;
    }

    // Frees the memory allocated for the architecture struct
    void FreeArchitecture(Architecture *architecture)
    {
        // freeing the cross links memory:
        for(int i = 0; i < architecture->numberOfCrossLinks; i++)
        {
            free(architecture->crossLinks[i]);
        }
        free(architecture->crossLinks);
    }

    // Methods that can be directly called by the user to initialize and fill the Architecture struct:

    // Initializes the Architecture struct and reads the diameter database
    // Stores the read diameter in the struct
    void InitializeArcAndReadDiameter(Architecture* archDiameter, int numberOfMonomers, char* architectureName)
    {
        // Initializing Cross-Links object to read dimaeters:
        InitializeArchitecture(archDiameter, architectureName);
        char* databaseFilePath;
        SetDiametersDatabaseFilePath(&databaseFilePath, numberOfMonomers);
        ReadDiameter(archDiameter, databaseFilePath);
        free(databaseFilePath);
    }

    // Initializes the Architecture struct and reads the cross-links database
    // Stores the read cross-links in the struct
    void InitializeArcAndReadCrossLinks(Architecture* archCrossLinks, int numberOfMonomers, char* architectureName)
    {
        InitializeArchitecture(archCrossLinks, architectureName);
        char* databaseFilePath;
        SetCrossLinksDatabaseFilePath(&databaseFilePath, numberOfMonomers);
        ReadCrossLinks(archCrossLinks, databaseFilePath);
        free(databaseFilePath);
    }

    // Initializes the Architecture struct and reads the diameter and axis lengths
    // Stores the read diameter and axis length values in the struct
    void InitializeArcConfinementDimensions(Architecture* archConfinement, int numberOfMonomers, char* architectureName)
    {
        InitializeArchitecture(archConfinement, architectureName);
        ReadDiameterAndAxisLength(archConfinement, numberOfMonomers);
    }

    // Initializes and fills the Architecture struct by reading cross links, diameter, and axis length from the databases
    void InitializeAndReadArc(Architecture* architecture, int numberOfMonomers, char* architectureName)
    {
        InitializeArchitecture(architecture, architectureName);
        char* databaseFilePath;
        SetCrossLinksDatabaseFilePath(&databaseFilePath, numberOfMonomers);
        ReadCrossLinks(architecture, databaseFilePath);
        free(databaseFilePath);
        SetDiametersDatabaseFilePath(&databaseFilePath, numberOfMonomers);
        ReadDiameter(architecture, databaseFilePath);
        free(databaseFilePath);
        SetAxisLengthDatabaseFilePath(&databaseFilePath, numberOfMonomers);
        ReadAxisLength(architecture, databaseFilePath);
        free(databaseFilePath);
    }

    // Returns whether the passed pair of monomer indices represent a pair of bonded monomers (within a polymer)
    // Accepts the cross links onject, number of monomers in a single polymer, and the two monomer indices (assuming counting starting from 1)
    bool AreMonomersBonded(Architecture archCrossLinks, int numberOfMonomers, int monomerIndex1, int monomerIndex2)
    {
        // First, checking if the monomers are from the same polymer:
        int polymerIndex1 = (monomerIndex1 - 1) / numberOfMonomers;
        int polymerIndex2 = (monomerIndex2 - 1) / numberOfMonomers;
        if(polymerIndex1 != polymerIndex1) // belong to different polymers
            return false;
        
        // Next, checking if the monomers are consecutive along contour
        int distance = abs(monomerIndex1 - monomerIndex2);
        if(distance == 1) // connectivity for linear polymer
            return true;
        // For a ring:
        if(strcmp(archCrossLinks.architectureName, "Arc_Lin")) // non-zero when polymer is not linear
        {
            distance = fmin(distance, numberOfMonomers - distance);
            if(distance == 1)
                return true;
        }

        // Finally checking for cross links:
        int reducedIndex1 = (monomerIndex1 - 1) % numberOfMonomers + 1;
        int reducedIndex2 = (monomerIndex2 - 1) % numberOfMonomers + 1;
        for(int i = 0; i < archCrossLinks.numberOfCrossLinks; i++)
        {
            if(reducedIndex1 == archCrossLinks.crossLinks[i][0] && reducedIndex2 == archCrossLinks.crossLinks[i][1])
                return true;
            // reversed order:
            if(reducedIndex2 == archCrossLinks.crossLinks[i][0] && reducedIndex1 == archCrossLinks.crossLinks[i][1])
                return true;
        }
        
        // Default case: none of the above conditions met:
        return false;

    }

    // Methods that treat the stretches between cross-linked monomers (both inclusive) as regions:

    // Counts the number of monomers between the cross-linked monomers (both inclusive)
    // Args:
    //      crossLink: an ordered pair (array) of integers that indicate the cross-linked pair of monomers
    //      numberOfMonomers: The total number of monomers in a single polymer
    // TODO: Merge this with regions.h methods
    int GetCrossLinkRegionCount(int* crossLink, int numberOfMonomers)
    {
        return modulo(crossLink[1]- crossLink[0], numberOfMonomers) + 1;
        // The modulo makes it so that this works even if the upper and lower indices are separated by the ori discontinuity i.e the 200-1 bond.
    }

    // Returns the next monomer along the polymer topology bonded to the current monomer.
    // This is useful when one wants to iterate over the indices of a subset of monomers that might contain the ori-discontinuity (200-1 bond)
    // Args:
    //      currentIndex: The index of the current monomer (indexing starts from 1)
    //      numberOfMonomers: The total number of monomers in a single polymer
    // Returns:
    //      the index of the next monomer along the polymer contour (indexing starts from 1)
    // TODO: Merge this with regions.h methods
    int GetNextMonomer(int currentIndex, int numberOfMonomers)
    {
        if(currentIndex % numberOfMonomers == 0) // The last monomer encountered
        {
            int monomerIndex = currentIndex - 1; // the index of current monomer when indexing starts from 0
            int polymerIndex = monomerIndex / numberOfMonomers; // floor divide; calculating polymer index so that the next monomer can also be calculated to be of the same polymer
            int nextMonomerIndex = modulo(monomerIndex + 1, numberOfMonomers) + numberOfMonomers * polymerIndex;
            return nextMonomerIndex + 1; // when indexing starts from 1
        }
        else
            return currentIndex + 1;
    }
#endif