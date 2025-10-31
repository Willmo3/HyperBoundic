//
// Created by will on 10/31/25.
//

#ifndef PDENCLOSE_LWRFLUX_H
#define PDENCLOSE_LWRFLUX_H
#include "FluxFunction.hpp"
#include "domains/Numeric.hpp"

/**
 * Light water reactor flux.
 * Procedure derived from Phocus.
 * @tparam T Numeric type flux function operates over.
 */
template<typename T>
requires Numeric<T>
class LwrFlux final : public FluxFunction<T> {
public:
    T flux(T value) override {
        return value * (value - 1);
    }
    T derivative_flux(T value) override {
        return value * -2 + 1;
    }
};

#endif //PDENCLOSE_LWRFLUX_H