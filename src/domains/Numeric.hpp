//
// Created by will on 9/17/25.
//

#ifndef PDENCLOSE_NUMERIC_H
#define PDENCLOSE_NUMERIC_H
#include <cstdint>

/**
 * Numeric data suitable for approximation schemes.
 */
template<typename T>
concept Numeric = requires(T a, uint32_t power, double scalar, std::ostream& out)
{
    /*
     * Numeric-numeric operations
     */
    a + a;
    a - a;
    a / a;
    a * a;
    // Will require tanh later when affine domain supports.
    // a.tanh();
    a.pow(power);
    a.abs();

    // TODO: log
    // TODO: sqrt

    /*
     * Numeric-scalar operations
     */
    a * scalar;
    a + scalar;
    a - scalar;
    a / scalar;
    a < scalar;
    a <= scalar;
    a > scalar;
    a >= scalar;

    /*
     * Other utility operations
     */
    out << a;
};

#endif //PDENCLOSE_NUMERIC_H