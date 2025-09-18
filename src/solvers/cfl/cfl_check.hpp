//
// Created by will on 9/18/25.
//

#ifndef PDEAPPROX_CFL_CHECK_H
#define PDEAPPROX_CFL_CHECK_H
#include <cstdint>

#include "../../domains/Numeric.hpp"
#include "../flux/FluxFunction.hpp"

/*
 * For hyperbolic schemes, c_max is defined as one.
 */
const uint32_t c_max = 1;

template<typename T>
requires Numeric<T>
bool cfl_check(FluxFunction<T> *f, T mesh_point, double delta_t, double delta_x) {
    // TODO: this cfl impl seems very conservative in my testing -- ensure acceptable.
    // std::cout << f->derivative_flux(mesh_point) * delta_t / delta_x << std::endl;
    return f->derivative_flux(mesh_point) * delta_t / delta_x <= c_max;
}

#endif //PDEAPPROX_CFL_CHECK_H