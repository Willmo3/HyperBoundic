//
// Created by will on 9/17/25.
//

#ifndef PDEAPPROX_REAL_H
#define PDEAPPROX_REAL_H
#include <cstdint>
#include <iosfwd>

/**
 * Wrapper class for real numbers, compliant with template requirements of abstract domain.
 */
class Real {
public:
    Real(double value);
    ~Real();
    double value() const;

    /*
     * Operations
     */
    Real operator+(const Real& right) const;
    Real operator-(const Real& right) const;
    Real operator*(const Real& right) const;
    Real operator/(const Real& right) const;

    Real tanh() const;
    Real pow(uint32_t power) const;
private:
    double _value;
};

std::ostream& operator<<(std::ostream& os, Real &rhs);


#endif //PDEAPPROX_REAL_H