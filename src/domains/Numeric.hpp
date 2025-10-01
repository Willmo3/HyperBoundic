//
// Created by will on 9/17/25.
//

#ifndef PDEAPPROX_NUMERIC_H
#define PDEAPPROX_NUMERIC_H

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

#endif //PDEAPPROX_NUMERIC_H