//
// Created by will on 11/4/25.
//

#ifndef PDENCLOSE_MATCH_SOLVER_H
#define PDENCLOSE_MATCH_SOLVER_H
#include "DifferenceSolver.hpp"
#include "LaxFriedrichsSolver.hpp"
#include "LeapfrogSolver.hpp"
#include "domains/Numeric.hpp"

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

#endif //PDENCLOSE_MATCH_SOLVER_H