//
// Created by will on 9/23/25.
//

#ifndef WAFFINE_WAFFINEFORM_H
#define WAFFINE_WAFFINEFORM_H
#include <unordered_map>
#include <vector>

#include "noise_symbol.h"
#include "Winterval.hpp"

// Author: Will Morris
// Credit to https://github.com/ogay/libaffa for cpp implementation tips.
// Credit to https://github.com/vanweric/Affapy for approximation strategies.

/**
 * Library for high performance affine form computations.
 * Affine forms are a symbolic representation of a
 * Form: x_hat + sum(real coefficient * noise_symbol).
 */
class WaffineForm {

public:
    /*
     * Constructors
     */
    /**
     * @param center Real number center for affine form.
     * @param starting_coeffs error coefficients to prime the affine form with
     */
    WaffineForm(double center, const std::unordered_map<noise_symbol_t, double> &starting_coeffs);
    /**
     * @param interval Interval to construct center, error points from.
     */
    explicit WaffineForm(const Winterval &interval);

    /*
     * Accessors
     */
    std::string to_string() const;
    double center() const;
    double radius() const;
    /**
     *
     * @param symbol Symbol to get coefficient of.
     * @return the coefficient if noise symbol represented in this affine form, NaN otherwise.
     */
    double coeff_of(noise_symbol_t symbol) const;
    Winterval to_interval() const;

    /*
     * Operations
     */

    /*
     * Unary operations
     */
    WaffineForm operator-() const;

    WaffineForm abs() const;

    /*
     * Binary affine operations
     */
    /**
     * Additive union of error symbols.
     */
    WaffineForm operator+(const WaffineForm &other) const;
    /**
     * Subtractive union of error symbols.
     */
    WaffineForm operator-(const WaffineForm &other) const;
    /**
     * Outer product of real values with others' errors. Additional error symbol added to account for loss of precision.
     */
    WaffineForm operator*(const WaffineForm &right) const;
    /**
     * Product of this affine form and the inverse of the rhs.
     */
    WaffineForm operator/(const WaffineForm &right) const;

    /**
     * @return The affine form multiplied by itself power times.
     */
    WaffineForm pow(uint32_t power) const;

    /*
     * Scalar arithmetic operations
     */
    WaffineForm operator*(double other) const;
    WaffineForm operator+(double other) const;
    WaffineForm operator-(double other) const;
    WaffineForm operator/(double other) const;
    /*
     * Scalar comparison operations
     * Note: these will be compared with the interval form of the system, as this reduces to concrete domain.
     */
    bool operator<(double other) const;
    bool operator>(double other) const;
    bool operator<=(double other) const;
    bool operator>=(double other) const;
private:
    /*
     * Non-affine function approximator helpers.
     */

    /**
     * Construct an affine approximation of this form applied to some non-affine function.
     * Form: (x_hat = alpha * x0 + zeta) + sum over errors (alpha * x_i * epsilon_i) + delta * epsilon_k
     * Generalizes first-order Chebyshev form alpha_x + beta (although we don't always use Chebyshev)
     * Implementation derived from Affapy.
     *
     * @param alpha Scalar multiplier for each error term.
     * @param zeta Constant offset for center term
     * @param delta Coefficient for new error term.
     * @return The affine approximation.
     */
    WaffineForm approximate_affine_form(double alpha, double zeta, double delta) const;
    /**
     * @return A new affine form representing the inverse of this form.
     * Approximate with mini-range.
     */
    WaffineForm inv() const;

    /*
     * Assorted helpers.
     */

    /**
     * @return A deep copy of this affine form.
     */
    WaffineForm clone() const;

    /*
     * Fields
     */

    /**
     * Center point for affine shape.
     */
    double _center;
    /**
     * Map of noise symbols to their coefficients.
     * Note that values will be heap allocated, so this data structure has a fixed size.
     */
    std::unordered_map<noise_symbol_t, double> _coefficients;
};

// Using reference to remove redundant copy.
// TODO: reflect this change elsewhere.
std::ostream& operator<<(std::ostream& os, WaffineForm rhs);


#endif //WAFFINE_WAFFINEFORM_H
