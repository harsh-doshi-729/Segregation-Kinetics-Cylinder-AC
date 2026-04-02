// This header file contains a struct for a generalised particle whose pair correlation can be found with respect to similar neighbours
// Other functions are included to calculate the pair correlation function

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <string.h>

#include "../System_File_Paths/system_file_paths.h" // Include the system file paths header for file paths
#include MONOMER_DISTRIBUTION // Include the monomer distribution header for necessary functions

#ifndef PAIR_CORRELATION_H
    #define PAIR_CORRELATION_H // include guard

    // Defining a struct for a generalised particle:
    struct GeneralisedParticle
    {
        int uid; // Unique identifier for the particle; does not need to be updated for each timestep
        double coordinates[3]; // Cartesian coordinates of the particle in 3D space
        char* description; // A unique description of the particle, e.g., its type or name; useful for debugging
        // The description does not need to be updated for each timestep
    };

    typedef struct GeneralisedParticle GeneralisedParticle; // Typedef for easier usage

    // Initializes and allocates memory for a new GeneralisedParticle that is returned
    // The arguments are used to initialize the new GeneralisedParticle. Coordinates are set to (0, 0, 0)
    // Args:
    // - uid: Unique identifier for the particle
    // - description: A string describing the particle, e.g., its type or name; useful for debugging
    // Returns: A pointer to the newly created GeneralisedParticle
    GeneralisedParticle GetNewGenParticle(int uid, const char* description)
    {
        // Allocate memory for the description and copy it
        char* desc = (char*)malloc(strlen(description) + 1);
        if (!desc) 
        {
            printf("Memory allocation failed for description in GeneralizedParticle failed.\n");
            exit(1);
        }
        strcpy(desc, description);

        // Declaring and initializing the new particle:
        GeneralisedParticle newParticle = { uid, {0, 0, 0}, desc };
        return newParticle;
        // This way of returning a struct by value is not ideal for big and bulky structs.
        // The struct has to be copied before the control can be passed back to the caller's context
        // The better way is to pass the pointer of the struct and assign the elements that way. But this might not be possible for constant elements.
        // Neither way is possible for constant elements in the struct when one wants to fill an array with dynamically allocated memory with such structs.
    }

    // Updates the coordinates of the given GeneralisedParticle with new coordinates
    // Args:
    // - particle: Pointer to the GeneralisedParticle whose coordinates are to be updated
    // - newCoordinates: Array of 3 doubles representing the new Cartesian coordinates of the particle
    void UpdateParticleCoordinates(GeneralisedParticle* particle, double newCoordinates[3])
    {
        for(int i = 0; i < 3; i++)
        {
            particle->coordinates[i] = newCoordinates[i];
        }
    }

    // Calculates the distance between two GeneralisedParticles
    // Args:
    // - particle1: Pointer to the first GeneralisedParticle
    // - particle2: Pointer to the second GeneralisedParticle
    double CalculateGenParticleDistance(const GeneralisedParticle* particle1, const GeneralisedParticle* particle2)
    {
        double distance = 0;
        for(int i = 0; i < 3; i++)
            distance += pow(particle1->coordinates[i] - particle2->coordinates[i], 2);
        return sqrt(distance);
    }

    // Updates the distance histogram for a list of GeneralisedParticles
    // Args:
    // - histogram: Pointer to the monomer_distribution_data struct that holds the histogram data to be eventually used for the pair correlation function
    // - particles: Pointer to the array of GeneralisedParticles
    void UpdateDistanceHistogram(monomer_distribution_data* histogram, const GeneralisedParticle* particles)
    {
        for(int i = 0; i < histogram->numberOfMonomers; i++)
        {
            for(int j = i + 1; j < histogram->numberOfMonomers; j++)
            {
                double distance = CalculateGenParticleDistance(&particles[i], &particles[j]);
                IncrementMonomerCounts(histogram, distance);
            }
        }
    }

    // Calculates the pair correlation function from the histogram data
    // Args:
    // - histogram: Pointer to the monomer_distribution_data struct that holds the histogram data to be eventually used for the pair correlation function
    // - pairCorrelationFunction: A double array that has already been allocated with a size equal to the number of bins of the histogram
    // - systemVolume: The volume that any GeneralisedParticles can access in the system
    // - numberOfTimeSteps: The number of time steps of the simulation of the system over which the histogram is constructed
    void CalculatePairCorrelationFunction(const monomer_distribution_data* histogram, double* pairCorrelationFunction, double systemVolume, int numberOfTimeSteps)
    {
        double averageDensity = histogram->numberOfMonomers / systemVolume; // Average density of the system
        // The normalization factor without the radius^2 in the denominator; common for all radii
        double preNormalization = histogram->numberOfMonomers * 4 * M_PI * histogram->binWidth * averageDensity * numberOfTimeSteps / 2;
        // Dividing by 2 because we have avoided double counting the monomers while adding the particle distances to the histogram
        for(int i = 0; i < histogram->numberOfBins; i++)
        {
            double radius = i * histogram->binWidth; // The radius of the shell being considered for g(r); aka left bin edge
            double normalizationFactor = preNormalization;
            if(i != 0) // Not normalizing by zero radius; ignoring the radius term
                normalizationFactor = normalizationFactor * pow(radius, 2);
            pairCorrelationFunction[i] = histogram->monomerCounts[i] / normalizationFactor;
        }
    }

    // Frees the memory allocated for a GeneralisedParticle
    // Args:
    // - particle: Pointer to the GeneralisedParticle to be freed
    void FreeGenParticle(GeneralisedParticle* particle)
    {
        free((void *)particle->description); // Free the allocated memory for the description
    }

#endif

// Algorithm for how such a pair correlation may proceed:
// 1. Initialize the structs for file reading 
// 2. At each timestep of the simulation, do the following:
//    a. Read the positions of the particles from the file
//    b. Assign the position of the generalisation of particle based on the system of particles as a struct; do this for each generalised particle
//    c. Calculate the distance between each pair of generalised particles
// 3. Store the distances in a histogram or a similar data structure and calculate the pair correlation function by normalising
// 
// I still have to assign the coordinates of the generalised particles for each timestep. (I will need to do this anyway in a separate code)
// Is this worth it to do it within this generalised pair correlation function rather than as a separate code?
// One advantage of this generalised framework is that I can reuse the same code for different types of particles.