//
// Created by will on 10/31/25.
//

#ifndef PDENCLOSE_BUCKLEYLEVERETTFLUX_H
#define PDENCLOSE_BUCKLEYLEVERETTFLUX_H
#include "FluxFunction.hpp"
#include "domains/Numeric.hpp"

template<typename T>
requires Numeric<T>
class BuckleyLeverett final : public FluxFunction<T> {
public:
    /**
     * x^2 / ((x^2) + (1/4 (1 - x)^2)
     * @param value value to substitute in for x. We will derive from the underlying discretization, the S function in formal notation.
     * @return the result of invoking the flux function with value.
     */
    T flux(T value) override {
        // Using intermediate value to avoid introducing new noise symbols.
        auto squared = value.pow(2);
        return squared / (squared + (value * -1 + 1).pow(2) * 0.25);
    }

    /**
     * Numerator: -8x^2 - x
     * Denominator: (5x^2 -2x +1)^2 = 25x^4 - 20x^3 + 14x^2 - 4x + 1
     * @param value Value to substittue in for x.
     * @return the result of invoking the flux function with value.
     */
    T derivative_flux(T value) override {
        // Manually computing powers wrt each other to maintain noise symbols between different forms.
        // Additionally, the fast descent of affine squaring has minimal benefit at a low power such as four.

        auto squared = value * value;
        auto cubed = squared * value;
        auto fourth = cubed * value;

        auto numer = squared * -8 - value;
        auto denom = fourth * 25 - cubed * 20 + squared * 14 - value * 4 + 1;

        return numer / denom;
    }
};

#endif //PDENCLOSE_BUCKLEYLEVERETTFLUX_H