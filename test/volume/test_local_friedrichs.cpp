//
// Created by will on 11/19/25.
//

#include <gtest/gtest.h>

#include "domains/Real.hpp"
#include "flux/BurgersFlux.hpp"
#include "meshes/CflCheck.hpp"
#include "solvers/volume/LocalLaxFriedrichsSolver.hpp"

TEST(llf, real_llf) {
    auto discretization_size = 5;
    auto num_timesteps = 4;

    auto initial_conditions = std::vector<Real>(discretization_size);
    initial_conditions[0] = 1.39;
    initial_conditions[1] = 2.66;
    initial_conditions[2] = 2.84;
    initial_conditions[3] = 2.75;
    initial_conditions[4] = 1.21;

    auto delta_t = 0.01;
    auto width_values = std::vector<double>(discretization_size);
    width_values[0] = 1;
    width_values[1] = 1;
    width_values[2] = 1;
    width_values[3] = 1;
    width_values[4] = 1;

    auto solver = LocalLaxFriedrichsSolver<Real>();
    auto solution_matrix = solver.solve(initial_conditions, width_values, discretization_size, num_timesteps, delta_t, new BurgersFlux<Real>);

    ASSERT_TRUE(solver.cfl_check_mesh(solution_matrix, new BurgersFlux<Real>, delta_t, width_values));
    ASSERT_NEAR(solution_matrix.get(3, 0).value(), 14.063770, 1e-5);
    ASSERT_NEAR(solution_matrix.get(3, 1).value(), 10.148959, 1e-5);
    ASSERT_NEAR(solution_matrix.get(3, 2).value(), 2.661306, 1e-5);
    ASSERT_NEAR(solution_matrix.get(3, 3).value(), 3.411484, 1e-5);
    ASSERT_NEAR(solution_matrix.get(3, 4).value(), 17.701213, 1e-5);
}