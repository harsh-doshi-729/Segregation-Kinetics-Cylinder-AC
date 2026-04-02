#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <stdbool.h>
#include "../System_File_Paths/system_file_paths.h"
#include RANDOM

void PrintRandomSeeds(int numberofSeeds)
{
    int limit = 900 * (int) (pow(10, 6)); // 900 million; limit for Marsaglia RNG
    //printing valid random integers (less than 900 million) in array format: 
    printf("[");
    for(int i = 0; i < numberofSeeds; i++)
    {   
        int seed = rand();
        if(seed >= limit)
        {
            i--; // reroll
            continue;
        }
        //else:
        if(i == numberofSeeds - 1) // last seed
            printf("%i]\n", seed);
        else
            printf("%i ", seed);
    }
}

void CorrectRandomSeeds(int numberofSeeds)
{
    //prior generated random seeds:
    int seeds[] = {1804289383, 846930886, 1681692777, 1714636915, 1957747793, 424238335, 719885386, 1649760492, 596516649, 1189641421, 1025202362, 1350490027, 783368690, 1102520059, 2044897763, 1967513926, 1365180540, 1540383426, 304089172, 1303455736, 35005211, 521595368, 294702567, 1726956429, 336465782};
    bool replacedSeeds[numberofSeeds];
    int threshold = 900 * (int) (pow(10, 6)); // 900 million
    //correcting the invalid seeds (for Marsaglia random number generator) and replacing them:
    for(int i = 0; i < numberofSeeds; i++)
    {
        if(seeds[i] > threshold) // invalid
        {
            replacedSeeds[i] = true;
            seeds[i] = generateRandomNumber(0, threshold);
        }
        else // valid
            replacedSeeds[i] = false;
    }

    //printing new seeds:
    printf("New seeds:\n[");
    for(int i = 0; i < numberofSeeds; i++)
    {   if(i == numberofSeeds - 1) // last seed
            printf("%i]\n", seeds[i]);
        else
            printf("%i ", seeds[i]);
    }

    //printing the indices that were replaced:
    printf("Replaced indices:\n[");
    for(int i = 0; i < numberofSeeds; i++)
    {
        if(replacedSeeds[i])
            printf("%i ", i);
    }
    printf("]\n");
}

int main(void)
{
    int numberofSeeds = 50;
    int metaSeed = 39; // Arc0 : 2, Arc2: 3, Arc-Loop4: 4, Arc-Loop10: 5, Arc-Loop20: 6, Arc0_relax_cyl: 7, Arc0_relax_no_confinement: 8, two Arc0 rings fixing+bonding: 9, Arc2-2: 10, Arc3: 11, Arc4: 12, Arc5: 13, Arc6: 14, Arc7: 15, Arc8: 16, Arc9: 17, Arc10: 27, Arc11: 19, relaxation cylinder_0.1: 20, relaxation cylinder_0.06: 21, relaxation cylinder_r5: 22, relaxation cylinder_r10: 23, fixing-bonding mixed states: 24, fbsr: 25, ArcI-8: 26, Arc1: 18, Arc1_1: 28, single_state_low: 29, single_state_high: 30, ArcR8_1_9: 31, Arc_Loop_5: 32, Arc_Lin: 33, ArcR8_1_19: 34, ArcR8_3_7: 35, Arc_Lin_Loop_1: 36, Arc_Lin_Loop_2: 37, Arc_Lin_Loop_5: 38, Arc_Lin_Loop_10: 39

    srand(metaSeed);

    // CorrectRandomSeeds(numberofSeeds);
    PrintRandomSeeds(numberofSeeds);

}
