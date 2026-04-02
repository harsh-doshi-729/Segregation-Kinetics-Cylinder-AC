#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#ifndef RANDOM_H // include guard
    #define RANDOM_H

    #ifndef M_PI
        #define M_PI 3.14159265358979323846
    #endif

    #define TOLERANCE 0.0000000001

    const double TWO_PI = 2 * M_PI;

    double generateRandomNumber(double lowerLimit, double upperLimit)
    {
        //set seed:
        //srand();
        double randomNumber = (double)rand()/(double)RAND_MAX;
            //printf(" %If\n", randomNumber);
            //a general formula to generate random numbers from a given range using Standard Uniform random numbers:
            return randomNumber * (upperLimit - lowerLimit) + lowerLimit;
    }

    void print10RandomNumbersToFile(FILE *randomNumbersFile)
    {
        for(int i = 0; i < 10; i++)
        {
            //printing to file:
            fprintf(randomNumbersFile, "%f ",generateRandomNumber(0, 1));//numbers separated by a space
        }
    }

    double getSumOfRandomNumbers(int sampleSize, int lowerLimit, int upperLimit )
    {
        double sum = 0;
        for(int i = 0; i < sampleSize; i++)
        {
            sum += generateRandomNumber(lowerLimit, upperLimit);
        }
        return sum;
    }

    int generateRandomWalkStep()
    {
        int remainder = rand() % 2;
        if(0 == remainder)
            return -1;
        else
            return 1;
    }

    //generates a pair of normal random variables by using the Box-Muller transform
    double* BoxMullerTransform(/*DEPRECATED: double uniform1, double uniform2, */double mean, double standardDeviation)//return a pointer to array
    {
        double uniform1 = 0;
        double uniform2 = 0;
        //excluding the value of 0 for our random variables
        // while(uniform1 <= TOLERANCE)
            uniform1 = generateRandomNumber(0, 1);
        // while(uniform2 <= TOLERANCE)
            uniform2 = generateRandomNumber(0, 1);
        //Box-Muller transform: optimized
        // ni = sqrt(-2*ln(R^2))*ui/R; R^2 = u1^2 + u2^2
        double *gaussianPair = malloc(sizeof (double) * 2);
        // double normSquared = pow(uniform1, 2) + pow(uniform2, 2);
        //checking if the numbers are within the unit circle
        // if(normSquared > 1)
            //return BoxMullerTransform(mean, standardDeviation); // rerolling the numbers; will this be a problem? Yes, I got sigkilled
        
        // double norm = sqrt(normSquared);
        //DEPRECATED:
        double firstFactor = sqrt( -2 * log(uniform1));
        gaussianPair[0] = firstFactor * cos(TWO_PI * uniform2) * standardDeviation - mean;
        gaussianPair[1] = firstFactor * sin(TWO_PI * uniform2) * standardDeviation - mean;
        // double firstFactor = sqrt( -2 * log(normSquared));
        // gaussianPair[0] = firstFactor * uniform1/norm;
        // gaussianPair[1] = firstFactor * uniform2/norm;
        //debugging:
        // if(isnan(gaussianPair[0])||isnan(gaussianPair[1]))
        //     printf("ERROR: Nan Generated\nNorm Squared = %lf\n", normSquared);
        return gaussianPair;
    }
#endif