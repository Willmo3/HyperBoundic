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
 * @return An initialized finite difference solver for a given domain
 */
template<typename T>
requires Numeric<T>
DifferenceSolver<T> *match_difference(const std::string &name) {
    if (name == "lax_friedrichs") {
        return new LaxFriedrichsSolver<T>();
    }
    if (name == "leapfrog") {
        return new LeapfrogSolver<T>();
    }
    std::cerr << "Unsupported PDE solvers!" << std::endl;
    exit(EXIT_FAILURE);
}


#endif //PDENCLOSE_MATCH_NAMES_H