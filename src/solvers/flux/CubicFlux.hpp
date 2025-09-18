//
// Created by will on 9/18/25.
//

#ifndef PDEAPPROX_CUBICFLUX_H
#define PDEAPPROX_CUBICFLUX_H

#include "FluxFunction.hpp"

/**
 *
 * @tparam T Numeric type flux function operates over.
 */
template<typename T>
requires Numeric<T>
class CubicFlux final : public FluxFunction<T> {
public:
    T flux(T value) override {
        return value.pow(3);
    }
    T derivative_flux(T value) override {
        return value.pow(2) * 3;
    }
};
#endif //PDEAPPROX_CUBICFLUX_H