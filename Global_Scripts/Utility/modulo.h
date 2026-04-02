// This is header file to implement a real modulo function for integers that works with negative operands as well
#ifndef MODULO_H
#define MODULO_H // an include guard; prevents compilation error when this file is included multiple times

// Performs the modulo of operand with respect to divisor. i.e. operand % divisor
int modulo(int operand, int divisor)
{
    const int result = operand % divisor;
    return result >= 0 ? result : result + divisor;
}

#endif