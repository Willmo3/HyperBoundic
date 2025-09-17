//
// Created by will on 9/17/25.
//

#ifndef PDEAPPROX_NUMERIC_H
#define PDEAPPROX_NUMERIC_H

/**
 * Numeric data suitable for approximation schemes.
 */
template<typename T>
concept Numeric = requires(T a, uint32_t power, std::ostream& out)
{
    a + a;
    a - a;
    a / a;
    a * a;
    a.tanh();
    a.pow(power);
    out << a;
};

#endif //PDEAPPROX_NUMERIC_H