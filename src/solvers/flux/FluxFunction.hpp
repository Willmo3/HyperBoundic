//
// Created by will on 9/18/25.
//

#ifndef PDENCLOSE_FLUXFUNCTION_H
#define PDENCLOSE_FLUXFUNCTION_H

#include "domains/Numeric.hpp"

template<typename T>
requires Numeric<T>
class FluxFunction {
public:
    FluxFunction() = default;
    virtual ~FluxFunction() = default;

    virtual T flux(T value) = 0;
    virtual T derivative_flux(T value) = 0;
};

#endif //PDENCLOSE_FLUXFUNCTION_H