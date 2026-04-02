#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <limits.h>

#ifndef MONOMER_DISTRIBUTION_H // include guard
    #define MONOMER_DISTRIBUTION_H

    #define TOLERANCE 0.0000000001
    // This header file contains structs and methods that are useful for calculating the position distribution of monomers/atoms in a system

    // A struct to store the relevant details required to calculate the monomer distribution
    struct monomer_distribution_data {
        int numberOfMonomers; // The total number of monomers in the system throughout the period during which the distribution is being calculated
        double lowerBoundary; // The position of the lower boundary for the monomer position
        double upperBoundary; // The position of the upper boundary for the monomer position
        int numberOfBins;
        double binWidth; // in the same units as the axisLength
        int* monomerCounts; // An array of counters to store the number of monomers occurring in each bin; aka the frequency distribution
        long int numberOfDataPoints; // The number of positions to make the entire distribution
        
    };

    // typedef for the struct:
    typedef struct monomer_distribution_data monomer_distribution_data;

    // INITIALIZATON:

    // Verifies that the value of the bond width and number of bins is valid
    static void VerifyBins(monomer_distribution_data *data)
    {
        if(data->binWidth <= 0 || data->numberOfBins <= 0)
        {
            printf("The number of bins and the bin width have not been set correctly! Terminating.\n");
            exit(1);
        }
    }

    // Verifies whether the value of the set dimensions is valid
    static void VerifyDimensions(monomer_distribution_data *data)
    {
        if(data->lowerBoundary >= data->upperBoundary)
        {
            printf("The boundaries of the simulation box have not been set correctly! Terminating.\n");
            exit(1);
        }
    }

    // Returns the length of the box along which the distribution is being calculated
    double GetAxisLength(monomer_distribution_data *data)
    {
        VerifyDimensions(data);
        return data->upperBoundary - data->lowerBoundary;
    }

    // Calculates lower and upper boundaries of the box based on the axis length and box centre
    void CalculateBoundaries(monomer_distribution_data *data, double axisLength, bool isBoxCentredAtOrigin)
    {
            data->lowerBoundary = 0;
            data->upperBoundary = axisLength;
            if(isBoxCentredAtOrigin)
            {
                data->lowerBoundary -= axisLength/2;
                data->upperBoundary -= axisLength/2;
            }
    }

    // Initializes the monomer distribution data with some constants.
    void InitializeMonomerDistributionData(monomer_distribution_data *data, int numberOfMonomers, double lowerBoundary, double upperBoundary)
    {
        // Note: the number of bins and the bin width is initialized separately. This is done so that the user has a choice as to what to initialize.
        data->numberOfMonomers = numberOfMonomers;
        data->lowerBoundary = lowerBoundary;
        data->upperBoundary = upperBoundary;
        data->numberOfDataPoints = 0;
    }

    // Allocates memory to and initializes the monomerCounts array. The number of bins must have already been set
    static void AllocateMonomerCountsArray(monomer_distribution_data *data)
    {
        VerifyBins(data);
        data->monomerCounts = (int*) calloc(data->numberOfBins, sizeof(*(data->monomerCounts)));
        assert(data->monomerCounts);
    }

    // static void FreeMonomerCounts(monomer_distribution_data *data)
    // {
    //     free(data->monomerCounts);
    // }


    // Rounds off the value of the axis length to the nearest multiple of the bin width
    static void RoundOffAxisLength(monomer_distribution_data *data, double binWidth)
    {
        double threshold = 0.00001;
        double axisLength = GetAxisLength(data);
        double remainder = fmod(axisLength, binWidth);
        if(remainder > threshold)
            axisLength = axisLength - remainder + binWidth; // rounding up true Axis Length to a multiple of binWidth
        // resetting the upper boundary:
        data->upperBoundary = data->lowerBoundary + axisLength;
        printf("Practical axis length has been set to %lf for bin width of %lf.\n", axisLength, binWidth);
    }

    // Accepts and sets the number of bins and sets the corresponding value of the bin width. The axis length must have already been set
    void SetNumberOfBins(monomer_distribution_data *data, int numberOfBins)
    {
        VerifyDimensions(data);
        data->numberOfBins = numberOfBins;
        data->binWidth = (GetAxisLength(data))/numberOfBins;
    }

    // Accepts and sets the bin width and sets the corresponding value of number of bins. The axis length must have already been set
    void SetBinWidth(monomer_distribution_data *data, double binWidth)
    {
        VerifyDimensions(data);
        RoundOffAxisLength(data, binWidth);
        data->binWidth = binWidth;
        data->numberOfBins = (GetAxisLength(data))/binWidth;
    }

    // Initializes the monomer_distribution data struct with the number of bins option
    void InitializeMonomerDistributionDataWithNumberOfBins(monomer_distribution_data *data, int numberOfMonomers, double lowerBoundary, double upperBoundary, int numberOfBins)
    {
        InitializeMonomerDistributionData(data, numberOfMonomers, lowerBoundary, upperBoundary);
        SetNumberOfBins(data, numberOfBins);
        AllocateMonomerCountsArray(data);
    }

    // Initializes the monomer_distribution data struct with the bin width option
    void InitializeMonomerDistributionDataWithBinWidth(monomer_distribution_data *data, int numberOfMonomers, double lowerBoundary, double upperBoundary, double binWidth)
    {
        InitializeMonomerDistributionData(data, numberOfMonomers, lowerBoundary, upperBoundary);
        SetBinWidth(data, binWidth);
        AllocateMonomerCountsArray(data);
    }

    // Initialize the monomer_distribution_data struct with the axis length and box centre
    void InitializeMonomerDistributionDataWithAxisLength(monomer_distribution_data *data, int numberOfMonomers, double axisLength, bool isBoxCentredAtOrigin)
    {
        data->numberOfMonomers = numberOfMonomers;
        CalculateBoundaries(data, axisLength, isBoxCentredAtOrigin);
        data->numberOfDataPoints = 0;
    }

    // Initialize the monomer distribution data struct with the axis length and number of bins
    void InitializeDistributionDataWithAxisLengthNumberOfBins(monomer_distribution_data *data, int numberOfMonomers, double axisLength, bool isBoxCentredAtOrigin, int numberOfBins)
    {
        InitializeMonomerDistributionDataWithAxisLength(data, numberOfMonomers, axisLength, isBoxCentredAtOrigin);
        SetNumberOfBins(data, numberOfBins);
        AllocateMonomerCountsArray(data);
    }

    // Initialize the monomer distribution data struct with the axis length and bin width
    void InitializeDistributionDataWithAxisLengthBinWidth(monomer_distribution_data *data, int numberOfMonomers, double axisLength, bool isBoxCentredAtOrigin, double binWidth)
    {
        InitializeMonomerDistributionDataWithAxisLength(data, numberOfMonomers, axisLength, isBoxCentredAtOrigin);
        SetBinWidth(data, binWidth);
        AllocateMonomerCountsArray(data);
    }

    // COMPUTATION:

    // Resets the values of the monomercounts array
    void ResetMonomerCounts(monomer_distribution_data *data)
    {
        for(int i = 0; i < data->numberOfBins; i++)
            data->monomerCounts[i] = 0;
    }

    // Updates the monomerCounts freqeuncy distribution with the position of a single monomer:
    void IncrementMonomerCounts(monomer_distribution_data *data, double position)
    {
        // Verifying if the position lies in range:
        if(position < data->lowerBoundary || position > data->upperBoundary)
        {
            // printf("The passed position is out of bounds! Terminating.\n");
            // exit(1);
            // Skipping the out of bounds data:
            return;
        }
        
        int index = floor((position - data->lowerBoundary)/data->binWidth);
        data->monomerCounts[index]++;
        data->numberOfDataPoints++;

        // Issuing a warning if the possibility of a count surpassing the capcity of an int variable is encountered; aka integer overflow:
        if(data->numberOfDataPoints >= INT_MAX)
        {
            printf("Warning: Due to the large number of data points, it is possible that integer overflow can occur for some of the bins in the distribution.\n");
        }
    }

    // Updates the monomerCounts frequency distribution with an array of positions (possibly the positions of an entire snapshot of the simulation):
    void UpdateMonomerCounts(monomer_distribution_data *data, double* positions, size_t arraySize)
    {
        for(int i = 0; i < arraySize; i++)
            IncrementMonomerCounts(data, positions[i]);
    }

    // Computes the final normalized distribution with a given normalization factor.
    // This can be useful if one wants to compute several distributions with a common normalization factor.
    void CalculateDistributionWithNormalizationFactor(monomer_distribution_data *data, double* distribution, size_t arrayLength, int normalizationFactor)
    {
        // matching the length of the array passed:
        if(arrayLength != data->numberOfBins)
        {
            printf("The length of the distribution array passed does not match the length of the frequency distribution! Terminating.\n");
            exit(1);
        }

        if(normalizationFactor == 0)
        {
            printf("The normalization factor passed to calculate the monomer distribution is 0! Divide by 0 error will occur. Terminating.\n");
            exit(1);
        }

        for(int i = 0; i < arrayLength; i++)
        {
            distribution[i] = (double)(data->monomerCounts[i]) / normalizationFactor / (data->binWidth); // normalizing with respect to the normalization factor and the bin width; This gives the probability distribution
            // Checking if the calculated distribution value is valid:
            if(distribution[i] < 0)
            {
                printf("The value of the distribution at index %i is negative (%lf)! Terminating.\n", i, distribution[i]);
                exit(1);
            }
        }
    }

    // Computes the final normalized distribution once all the positions have been added to the frequency distribution
    void CalculateDistribution(monomer_distribution_data *data, double* distribution, size_t arrayLength)
    {
        if(data->numberOfDataPoints == 0)
        {
            printf("No data points were added to the distribution! The distribution is empty. Terminating.\n");
            exit(1);
        }
        CalculateDistributionWithNormalizationFactor(data, distribution, arrayLength, data->numberOfDataPoints); // Normalizing with respect to the number of data points
        // This should give the ideal probability distribution
    }


    // Prints the distribution passed to the passed file pointer
    // If truncateDistribution is true, the 0 valued bins at the ends are skipped (if any are present)
    void PrintDistribution(monomer_distribution_data *data, double* distribution, bool truncateDistribution, FILE* writeFilePointer)
    {
        // TODO: figure out the ideal way to pass such arrays: pass the variable length array as it is

        bool finiteValueFound = false; //  A flag to indicate that the first finite distribution value from the left has been found
        bool rightZeroFound = false; // A flag to indicate that the first zero on the right side of the finite distribution has been found
        int leftZeroIndex = 0;
        int rightZeroIndex = data->numberOfBins - 1; // The indices of the zeroes flanking either side of the distribution.
        // The data outside this index range just contains zeroes, so it can be ignored

        if(truncateDistribution)
        {
            for(int i = 0; i < data->numberOfBins; i++) // Finding left zero index
            {
                if(distribution[i] > TOLERANCE)
                {
                    if(i == 0) // first element on the left is finite
                        leftZeroIndex = i;
                    else
                        leftZeroIndex = i - 1;
                    break;
                }

            }

            for(int i = data->numberOfBins - 1; i >= 0; i--) // finding right zero index
            {
                if(distribution[i] > TOLERANCE)
                {
                    if(i == data->numberOfBins-1) // last element on right side is finite
                        rightZeroIndex = i;
                    else
                        rightZeroIndex = i + 1;
                    break;
                }
            }
        }
        // Printing the distribution:
        fprintf(writeFilePointer, "Left Bin Edge, Probability density\n"); // The header
        for(int i = leftZeroIndex; i <= rightZeroIndex; i++)
        {
            double leftBinEdge = i * data->binWidth + data->lowerBoundary;
            fprintf(writeFilePointer, "%.6lf, %.6lf\n", leftBinEdge, distribution[i]);
        }
    }

    // Checks whether two double values are close enough to be considered equal
    bool AreAlmostEqual(double x, double y)
    {
        // taken from PurpleOrangeSkies' answer on this reddit thread: https://www.reddit.com/r/C_Programming/comments/4thsn7/comparing_doubles_in_c/?rdt=46060
        return (fabs(x - y) < (__DBL_EPSILON__));
    }

    // Checks whether two monomer distrbution objects are compatible for addition of monomer counts
    // Checks equality of the following variables: numberOfBins, binWdith, lowerBoundary, upperBoundary
    bool AreMonomerDistributionsSummable(monomer_distribution_data dist_data1, monomer_distribution_data dist_data2)
    {
        bool areSummable = true;
        areSummable *= (dist_data1.numberOfBins == dist_data2.numberOfBins);
        areSummable *= AreAlmostEqual(dist_data1.binWidth, dist_data2.binWidth);
        areSummable *= AreAlmostEqual(dist_data1.lowerBoundary, dist_data2.lowerBoundary);
        areSummable *= AreAlmostEqual(dist_data1.upperBoundary, dist_data2.upperBoundary);
        return areSummable;
    }

    // Checks whether the properties (member variables) of the two passed monomer dist objects are the same.
    // Checks equality of the following variables: numberOfMonomers, numberOfBins, binWdith, lowerBoundary, upperBoundary
    bool AreMonomerDistributionsEqual(monomer_distribution_data dist_data1, monomer_distribution_data dist_data2)
    {
        return (dist_data1.numberOfMonomers == dist_data2.numberOfMonomers) && AreMonomerDistributionsSummable(dist_data1, dist_data2);
    }


    // Adds the monomer counts of two different passed monomder_distribution_data objects and stores it into a third (uninitialized) object
    // Only adds the two if the bin width, number of bins, lower and upper boundaries match for both objects
    // It assumes the data used for the two distributions are mutually exclusive
    void AddMonomerCounts(monomer_distribution_data dist_data1, monomer_distribution_data dist_data2, monomer_distribution_data* result_dist_data)
    {
        // Verifying that the number of bins and bin width is same for both objects:
        if(!AreMonomerDistributionsSummable(dist_data1, dist_data2))
        {
            printf("ERROR: Adding monomer counts of the two distribution objects failed since the they were not summable!\n");
            exit(1);
        }

        // Once verified:
        // Initializing the resultant object to match the sum properties
        int totalMonomers =  dist_data1.numberOfMonomers + dist_data2.numberOfMonomers;
        InitializeMonomerDistributionDataWithBinWidth(result_dist_data, totalMonomers, dist_data1.lowerBoundary, dist_data1.upperBoundary, dist_data1.binWidth);
        ResetMonomerCounts(result_dist_data);
        for(int i = 0; i < result_dist_data->numberOfBins; i++)
        {
            result_dist_data->monomerCounts[i] = dist_data1.monomerCounts[i] + dist_data2.monomerCounts[i];
        }
        result_dist_data->numberOfDataPoints = dist_data1.numberOfDataPoints + dist_data2.numberOfDataPoints;
    }

    // Frees the allocated memory for calculating the frequency distribution:
    void FreeMonomerDistribution(monomer_distribution_data *data)
    {
        free(data->monomerCounts);
    }
#endif