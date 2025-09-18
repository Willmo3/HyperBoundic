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
     * NOTE: do not use reference so that temporary values can be applied in the way e.g. doubles are.
     */
    Real operator+(Real right) const;
    Real operator-(Real right) const;
    Real operator*(Real right) const;
    Real operator/(Real right) const;
    bool operator==(Real right) const;

    Real tanh() const;
    Real pow(uint32_t power) const;
private:
    double _value;
};

// Note: cannot use reference for rhs because we want to be able to print shortlived values
// i.e. std::cout << Real(a) + Real(b) << std::endl;
std::ostream& operator<<(std::ostream& os, Real rhs);


#endif //PDEAPPROX_REAL_H