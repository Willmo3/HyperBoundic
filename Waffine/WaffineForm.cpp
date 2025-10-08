//
// Created by will on 9/23/25.
//

#include "WaffineForm.hpp"

#include <cmath>
#include <cstring>
#include <format>
#include <iostream>
#include <map>
#include <numeric>
#include <ranges>
#include <vector>

/*
 * Constructors
 */
WaffineForm::WaffineForm(double center, const std::unordered_map<noise_symbol_t, double> &starting_coeffs):
    _center(center), _coefficients(std::unordered_map<noise_symbol_t, double>()) {
    // Initialize this map with the explicitly defined starting values.
    for (auto pair : starting_coeffs) {
        _coefficients.insert(pair);
    }
}

WaffineForm::WaffineForm(const Winterval& interval): _center((interval.min() + interval.max()) / 2),
    _coefficients(std::unordered_map<noise_symbol_t, double>()) {
    _coefficients.insert(std::pair(new_noise_symbol(), (interval.min() - interval.max()) / 2));
}

/*
 * Unary operators
 */
WaffineForm WaffineForm::operator-() const {
    auto value = clone();
    value._center = -_center;
    for (const auto symbol: _coefficients | std::views::keys) {
        value._coefficients[symbol] *= -1;
    }
    return value;
}
WaffineForm WaffineForm::abs() const {
    // Strictly negative
    if (this->operator<(0)) {
        return this->operator*(-1);
    }
    // Strictly positive
    if (this->operator>(0)) {
        return clone();
    }
    // Straddles
    // Return abs(this.center / 2) + sum (this.noise / 2)
    auto value = clone();
    value._center = std::abs(value._center / 2);
    for (auto symbol: value._coefficients | std::views::keys) {
        value._coefficients[symbol] /= 2;
    }
    return value;
}

/*
 * Affine arithmetic operators.
 */
WaffineForm WaffineForm::operator+(const WaffineForm &other) const {
    auto value = clone();
    value._center += other._center;

    for (auto [symbol, coeff] : other._coefficients) {
        // Outer product is union of both fields' error symbols. Common error symbols are added.
        if (!value._coefficients.contains(symbol)) {
            value._coefficients[symbol] = coeff;
        } else {
            value._coefficients[symbol] += coeff;
        }
    }
    // Since affine addition introduces no new error, we don't need to add a new value!
    return value;
}
WaffineForm WaffineForm::operator-(const WaffineForm &other) const {
    auto value = clone();
    value._center -= other._center;

    for (auto [symbol, coeff] : other._coefficients) {
        if (!value._coefficients.contains(symbol)) {
            value._coefficients[symbol] = -coeff;
        } else {
            value._coefficients[symbol] -= coeff;
        }
    }

    return value;
}
WaffineForm WaffineForm::operator*(const WaffineForm &right) const {
    auto result = WaffineForm(this->_center * right._center, std::unordered_map<noise_symbol_t, double>());
    // Affine form multiplication is an outer product.

    // Perform product for all error symbols in rhs.
    for (auto [symbol, coeff] : right._coefficients) {
        if (!this->_coefficients.contains(symbol)) {
            // Add missing error terms scaled by left's center.
            result._coefficients[symbol] = this->_center * coeff;
        } else if (this->_coefficients.contains(symbol)) {
            auto a = this->_center * coeff;
            auto b = right._center * this->_coefficients.at(symbol);
            result._coefficients[symbol] = a + b;
        }
    }

    // Now go and perform similar calculation for error symbols in lhs that weren't caught earlier.
    for (auto [symbol, coeff] : this->_coefficients) {
        if (!right._coefficients.contains(symbol)) {
            result._coefficients[symbol] = right._center * coeff;
        }
    }

    // Affine multiplication adds a noise symbol.
    // For now, we add an error w/ coeff rad * rad, following Affapy impl.
    result._coefficients[new_noise_symbol()] = this->radius() * right.radius();
    return result;
}
WaffineForm WaffineForm::operator/(const WaffineForm &right) const {
    return operator*(right.inv());
}

/*
 * Scalar arithmetic operators
 */
WaffineForm WaffineForm::operator*(double other) const {
    auto value = clone();
    value._center *= other;
    for (const auto symbol: _coefficients | std::views::keys) {
        value._coefficients[symbol] *= other;
    }
    return value;
}
WaffineForm WaffineForm::operator+(double other) const {
    auto value = clone();
    value._center += other;
    // Notice: addition does not affect error symbols.
    // effectively, the polytope is simply being translated.
    return value;
}
WaffineForm WaffineForm::operator-(double other) const {
    auto value = clone();
    value._center -= other;
    return value;
}
WaffineForm WaffineForm::operator/(double other) const {
    // Special case: 0. In affine forms, this will set all terms to 0, leading to a unit form.
    if (other == 0) {
        return { 0, std::unordered_map<noise_symbol_t, double>() };
    }
    return operator*(1 / other);
}
WaffineForm WaffineForm::pow(uint32_t power) const {
    // TODO: as we adapt the numeric API, we could switch this to use negative numbers w/ the inverse strategy.
    if (power == 0) {
        // Our implementation always returns affine forms, even if the power is 0 -- in this case, an exact affine form.
        return { 1, std::unordered_map<noise_symbol_t, double>() };
    }

    auto odd_power = power > 1 && power % 2 == 1;
    if (odd_power) {
        // Can descend logarithmically given an even power. We will do the extra mult later.
        power -= 1;
    }

    auto result = clone();
    while (power > 1) {
        // Perform multiply and half power each time until down to pow 1, unit operation.
        // Insight: squaring intermediate results allows our quick logarithmic descent.
        result = result * result;
        power /= 2;
    }

    if (odd_power) {
        // Now perform that last standard multiplication we saved.
        return result * *this;
    }
    return result;
}

/*
 * Scalar comparison operators.
 */
bool WaffineForm::operator<(double other) const {
    return to_interval() < other;
}
bool WaffineForm::operator>(double other) const {
    return to_interval() > other;
}
bool WaffineForm::operator<=(double other) const {
    return to_interval() <= other;
}
bool WaffineForm::operator>=(double other) const {
    return to_interval() >= other;
}

/*
 * Accessors
 */
std::string WaffineForm::to_string() const {
    std::string retval = std::string();
    retval += "Interval concretization: ";
    retval += "[" + std::to_string(to_interval().min());
    retval += ", ";
    retval += std::to_string(to_interval().max()) + "]\n";
    retval += "Center: " + std::to_string(_center) + "\n";
    retval += "Radius: " + std::to_string(radius()) + "\n";
    retval += "Noise symbols:";
    for (auto [symbol, coeff] : _coefficients) {
        retval += " (" + std::to_string(symbol) + ": " + std::to_string(coeff) + "),";
    }
    retval += "\n";
    return retval;
}

double WaffineForm::center() const {
    return _center;
}

double WaffineForm::radius() const {
    return std::accumulate(_coefficients.begin(), _coefficients.end(), 0.0,
        [](auto sum, auto pair) { return sum + std::abs(pair.second); });
}

double WaffineForm::coeff_of(noise_symbol_t symbol) const {
    if (_coefficients.contains(symbol)) {
        return _coefficients.at(symbol);
    }
    return NAN;
}

Winterval WaffineForm::to_interval() const {
    // Note: unable to do accumulation because of behavior with unordered maps.
    double error_magnitude = 0;
    for (auto coeff: _coefficients | std::views::values) {
        error_magnitude += std::abs(coeff);
    }
    return {_center - error_magnitude, _center + error_magnitude};
}

/*
 * Non-affine approximators
 */

/*

Approximating a non-affine form follows a general pattern:
- Create a new center with new center alpha * old_center + zeta
- Scale each error term by alpha
- Add a new error term with coeff delta.

*/
WaffineForm WaffineForm::approximate_affine_form(double alpha, double zeta, double delta) const {
    auto center = alpha * _center + zeta;

    auto map = std::unordered_map<noise_symbol_t, double>();
    for (auto [symbol, coeff] : _coefficients) {
        auto new_value = alpha * coeff;
        map.insert({symbol, new_value});
    }

    map.insert({new_noise_symbol(), delta});
    return { center, map };
}

/*
 * Use mini-range approximation rather than standard Chebyshev.
 * Credit: libaffa
 */
WaffineForm WaffineForm::inv() const {
    auto interval = this->to_interval();
    if (interval.contains(0)) {
        // If interval contains 0, infinity will be in this interval, or the interval will just be [0, 0].
        // The blowup to infinity wipes away the dependence of the variables and adds a term of unlimited magnitude.
        // We will hence construct a new affine form from the corresponding interval, since no dependencence can be preserved.
        return WaffineForm(this->to_interval());
    }

    auto a = interval.abs().min();
    auto b = interval.abs().max();
    // Derivative of 1/x: -1/x^2. Notice: series defined over first derivative.
    auto alpha = -1 / std::pow(b, 2);

    auto range = Winterval((1/a - alpha * a), 2 / b);
    auto zeta = range.mid();

    // if negative value included in interval, flip result.
    if (interval.min() < 0) {
        zeta = -zeta;
    }

    // New noise term will be radius of mini-range interval.
    return approximate_affine_form(alpha, zeta, range.radius());
}

/*
 * Internal helpers
 */
// Internal clone constructor
WaffineForm WaffineForm::clone() const {
    return { this->_center, std::unordered_map(this->_coefficients) };
}

/*
 * Associated operators.
 */
std::ostream& operator<<(std::ostream& os, WaffineForm rhs) {
    os << rhs.to_string();
    return os;
}