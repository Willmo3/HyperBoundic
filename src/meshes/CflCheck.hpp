//
// Created by will on 9/18/25.
//

#ifndef PDENCLOSE_CFL_CHECK_H
#define PDENCLOSE_CFL_CHECK_H
#include <cstdint>

#include "domains/Numeric.hpp"
#include "flux/FluxFunction.hpp"

/*
 * For explicit schemes, c_max is defined as one.
 */
const uint32_t c_max = 1;

template<typename T>
requires Numeric<T>
bool cfl_check(FluxFunction<T> *f, T mesh_point, double delta_t, double delta_x) {
    auto cfl_value = f->derivative_flux(mesh_point).abs() * delta_t / delta_x;
    if (cfl_value >= c_max) {
        std::cerr << "CFL check failed with value: " << cfl_value << std::endl;
        return false;
    }
    return true;
}

#endif //PDENCLOSE_CFL_CHECK_H