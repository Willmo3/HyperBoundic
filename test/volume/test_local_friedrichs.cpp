//
// Created by will on 11/19/25.
//

#include <gtest/gtest.h>

#include "domains/Real.hpp"
#include "flux/BurgersFlux.hpp"
#include "solvers/volume/LocalLaxFriedrichsSolver.hpp"

TEST(llf, real_llf) {
    auto discretization_size = 4;
    auto num_timesteps = 4;

    auto initial_conditions = std::vector<Real>(discretization_size);
    initial_conditions[0] = 1.0;
    initial_conditions[1] = 2.0;
    initial_conditions[2] = 3.0;
    initial_conditions[3] = 4.0;

    auto delta_t = 0.02;
    auto width_values = std::vector<double>(discretization_size);
    width_values[0] = 100;
    width_values[1] = 100;
    width_values[2] = 100;
    width_values[3] = 100;

    auto solution_matrix = LocalLaxFriedrichsSolver<Real>().solve(initial_conditions, width_values, discretization_size, num_timesteps, delta_t, new BurgersFlux<Real>);
    // solution_matrix.print_system();
}