//
// Created by will on 11/4/25.
//

#ifndef PDENCLOSE_SELECT_FLUX_H
#define PDENCLOSE_SELECT_FLUX_H
#include <iostream>
#include <string>

#include "BuckleyLeverettFlux.hpp"
#include "BurgersFlux.hpp"
#include "CubicFlux.hpp"
#include "FluxFunction.hpp"
#include "LwrFlux.hpp"
#include "domains/Numeric.hpp"

/**
 *
 * @tparam T Numeric type to operate over.
 * @param name Name of the flux function
 * @return An initialized flux function corresponding to the name.
 */
template<typename T>
requires Numeric<T>
FluxFunction<T> *match_flux(const std::string &name) {
    if (name == "buckley_leverett") {
        return new BuckleyLeverett<T>();
    }
    if (name == "burgers") {
        return new BurgersFlux<T>();
    }
    if (name == "cubic") {
        return new CubicFlux<T>();
    }
    if (name == "lwr") {
        return new LwrFlux<T>();
    }
    std::cerr << "Unsupported flux function!" << std::endl;
    exit(EXIT_FAILURE);
}

#endif //PDENCLOSE_SELECT_FLUX_H
