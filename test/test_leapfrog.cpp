//
// Created by will on 10/17/25.
//

#include "../lib/Waffine/WaffineForm.hpp"
#include "../src/domains/Real.hpp"
#include "../src/solvers/LeapfrogSolver.hpp"
#include "../src/solvers/flux/CubicFlux.hpp"
#include "../lib/Winterval/Winterval.hpp"
#include "../lib/Wixed/WixedForm.hpp"
#include "gtest/gtest.h"

TEST(leapfrog, real_approx) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_t = 0.02; // Spacing of time. Assuming operating over 4 logical time split into 4 timesteps = 1.
    double delta_x = 1; // Spatial discretization. Assuming total space of four split into 4 parts = 1.

    auto initial_conditions = std::unique_ptr<Real>(static_cast<Real *>(calloc(discretization_size, sizeof(double))));

    initial_conditions.get()[0] = 1.0;
    initial_conditions.get()[1] = 2.0;
    initial_conditions.get()[2] = 3.0;
    initial_conditions.get()[3] = 4.0;

    auto solution_matrix = LeapfrogSolver<Real>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<Real>());
    // test the third row of the discretization.
    ASSERT_NEAR(solution_matrix.get(2, 0).value(), 1.125503, 0.001);
    ASSERT_NEAR(solution_matrix.get(2, 1).value(), 2.611825, 0.001);
    ASSERT_NEAR(solution_matrix.get(2, 2).value(), 2.874497, 0.001);
    ASSERT_NEAR(solution_matrix.get(2, 3).value(), 3.388175, 0.001);
}
