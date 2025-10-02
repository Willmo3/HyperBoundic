//
// Created by will on 9/17/25.
//

#include <cmath>

#include "Real.hpp"

#include <iosfwd>
#include <string>

Real::Real(double value): _value(value) {}
Real::~Real() = default;

double Real::value() const {
    return _value;
}

/*
 * Operations
 */
Real Real::operator+(const Real &right) const {
    return { _value + right._value };
}
Real Real::operator-(const Real &right) const {
    return { _value - right._value };
}
Real Real::operator*(const Real &right) const {
    return { _value * right._value };
}
Real Real::operator/(const Real &right) const {
    return { _value / right._value };
}
bool Real::operator==(const Real &right) const {
    return _value == right._value;
}
bool Real::operator<=(const Real &right) const {
    return _value <= right._value;
}
bool Real::operator<(const Real &right) const {
    return _value < right._value;
}
bool Real::operator>=(const Real &right) const {
    return _value >= right._value;
}
bool Real::operator>(const Real &right) const {
    return _value > right._value;
}

std::ostream& operator<<(std::ostream& os, Real rhs) {
    os << std::to_string(rhs.value());
    return os;
}

Real Real::tanh() const {
    return { std::tanh(_value) };
}
Real Real::pow(uint32_t power) const {
    return { std::pow(_value, power) };
}
Real Real::abs() const {
    return { std::abs(_value) };
}
