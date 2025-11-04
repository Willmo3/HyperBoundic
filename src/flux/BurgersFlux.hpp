//
// Created by will on 10/31/25.
//

#ifndef PDENCLOSE_BURGERSFLUX_H
#define PDENCLOSE_BURGERSFLUX_H
#include "FluxFunction.hpp"
#include "domains/Numeric.hpp"

template<typename T>
requires Numeric<T>
class BurgersFlux final : public FluxFunction<T> {
public:
    T flux(T value) override {
        return value.pow(2) * 0.5;
    }
    T derivative_flux(T value) override {
        return value;
    }
};

#endif //PDENCLOSE_BURGERSFLUX_H