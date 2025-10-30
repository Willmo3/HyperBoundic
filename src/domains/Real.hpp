//
// Created by will on 9/17/25.
//

#ifndef PDENCLOSE_REAL_H
#define PDENCLOSE_REAL_H
#include <cstdint>
#include <iosfwd>

#include "cereal/archives/binary.hpp"


/**
 * Wrapper class for real numbers, compliant with template requirements of abstract domain.
 */
class Real {
public:
    Real() = default;
    Real(double value);
    ~Real();
    double value() const;

    /*
     * Operations
     */
    Real operator+(const Real &right) const;
    Real operator-(const Real &right) const;
    Real operator*(const Real &right) const;
    Real operator/(const Real &right) const;
    bool operator==(const Real &right) const;
    bool operator<=(const Real &right) const;
    bool operator<(const Real &right) const;
    bool operator>=(const Real &right) const;
    bool operator>(const Real &right) const;

    Real tanh() const;
    Real pow(uint32_t power) const;
    Real abs() const;

    /*
     * Serialization support through cereal.
     */
    template<class Archive>
    void serialize(Archive & archive) {
        archive(_value);
    }
private:
    double _value;
};

// Note: cannot use reference for rhs because we want to be able to print shortlived values
// i.e. std::cout << Real(a) + Real(b) << std::endl;
std::ostream& operator<<(std::ostream& os, Real rhs);


#endif //PDENCLOSE_REAL_H