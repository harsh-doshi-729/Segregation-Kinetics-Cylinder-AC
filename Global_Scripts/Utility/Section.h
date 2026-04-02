// This header file contains variables and functions to divide the simulation box into multiple sections
// The aim is to perform analysis on these sections separately to compare across sections

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#ifndef SECTION_H // include guard
    #define SECTION_H
    const int numberOfSections = 8;

    // Returns the section that the passed monomer belongs to based on its position, the number of monomers in a single polymer, radius of the confinement (without wall), and total no. of sections
    int GetSectionID(int numberOfMonomers, double radius, int numberOfSections, double* position)
    {
        // The division of the cylinder into different sections will be relevant here. 
        // As of now, the sections are made in the x-y plane and projected along the z axis.
        // The monomer distribution is calculated along the z-axis
        // All these sections must have the same volume
        
        // The x-y circular cross-section is divided into 4 sections of a circle.
        // Further, each section is divided along the radial direction into a smaller section and section of an annular ring
        // The smaller sections are indexed first with increasing theta (anticlockwise) starting from the -ve x-axis. Then the outer sections 
        // are indexed in a similar fashion
        double wall = 0.5; // The thickness of the wall in sigma = 1 units
        int numberOfAngularSections = numberOfSections / 2;
        double radialThreshold = (radius - wall)/ sqrt(2); // The threshold between the smaller section and the outer sections

        // calculating radial and theta coordinate of the position:
        double r = sqrt(pow(position[0], 2) + pow(position[1], 2));
        if(position[0] == 0)
            position[0] += 0.000001; // handling zero case specially
        double theta = atan2(position[1], position[0]); // arctan(y/x) in radian from -Pi to +Pi
        theta += M_PI; // from 0 to 2*PI
        double sectionWidth = 2 * M_PI / numberOfAngularSections; // Angular width in radian
        int sectionID = floor(theta / sectionWidth) + 1; // The indexing starts from 1
        if(r > radialThreshold)
            sectionID += numberOfAngularSections;

        // Debugging:
        // printf("r: %lf\ntheta: %lf\n", r, theta);

        return sectionID;
    }
#endif