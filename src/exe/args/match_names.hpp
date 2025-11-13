//
// Created by will on 11/4/25.
//

#ifndef PDENCLOSE_MATCH_NAMES_H
#define PDENCLOSE_MATCH_NAMES_H
#include <cstdlib>
#include <string>

#include "flux/BuckleyLeverettFlux.hpp"
#include "flux/BurgersFlux.hpp"
#include "flux/CubicFlux.hpp"
#include "flux/FluxFunction.hpp"
#include "flux/LwrFlux.hpp"
#include "solvers/difference/LeapfrogSolver.hpp"

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

/**
 *
 * @tparam T Numeric type to operate over
 * @param name Name of finite difference scheme to match
 * @param initial_state Initial state of the system
 * @param discretization_size Size of spatial discretization
 * @param num_timesteps Number of timesteps for the approximation
 * @param delta_t Time discretization
 * @param delta_x Space discretization
 * @param flux Flux function to use for this solver
 * @return An initialized finite difference solver for a given domain
 */
template<typename T>
requires Numeric<T>
// TODO: add new fields to generate solvers
DifferenceSolver<T> *match_difference(const std::string &name,
    const std::vector<T> &initial_state,
    uint32_t discretization_size,
    uint32_t num_timesteps,
    double delta_t,
    double delta_x,
    FluxFunction<T> *flux) {

    if (name == "lax_friedrichs") {
        return new LaxFriedrichsSolver<T>(initial_state, discretization_size, num_timesteps, delta_t, delta_x, flux);
    }
    if (name == "leapfrog") {
        return new LeapfrogSolver<T>(initial_state, discretization_size, num_timesteps, delta_t, delta_x, flux);
    }
    std::cerr << "Unsupported PDE solvers!" << std::endl;
    exit(EXIT_FAILURE);
}


#endif //PDENCLOSE_MATCH_NAMES_H