//
// Created by will on 9/18/25.
//

#include <gtest/gtest.h>
#include "../src/domains/Real.hpp"
#include "../src/solvers/LaxFriedrichsSolver.hpp"
#include "../src/solvers/flux/CubicFlux.hpp"
#include "../Winterval/src/Winterval.h"

TEST(friedrichs, real_approx) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_x = 1; // Spatial discretization. Assuming total space of four split into 4 parts = 1.
    double delta_t = 1; // Spacing of time. Assuming operating over 4 logical time split into 4 timesteps = 1.

    auto initial_conditions = static_cast<Real *>(calloc(discretization_size, sizeof(double)));
    assert(initial_conditions);

    initial_conditions[0] = 1.0;
    initial_conditions[1] = 2.0;
    initial_conditions[2] = 3.0;
    initial_conditions[3] = 4.0;

    auto solution_matrix = LaxFriedrichsSolver<Real>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<Real>());
    // test the third row of the discretization.
    ASSERT_EQ(solution_matrix.get(2, 0), 2355);
    ASSERT_EQ(solution_matrix.get(2, 1), 22711);
    ASSERT_EQ(solution_matrix.get(2, 2), -2351);
    ASSERT_EQ(solution_matrix.get(2, 3), -22705);

    free(initial_conditions);
    initial_conditions = nullptr;
}

TEST(friedrichs, interval_approx) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_x = 1;
    double delta_t = 1;

    auto initial_conditions = static_cast<Winterval *>(calloc(discretization_size, sizeof(Winterval)));
    assert(initial_conditions);

    initial_conditions[0] = Winterval(0, 1);
    initial_conditions[1] = Winterval(1, 2);
    initial_conditions[2] = Winterval(2, 3);
    initial_conditions[3] = Winterval(3, 4);

    auto solution_matrix = LaxFriedrichsSolver<Winterval>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<Winterval>());
    ASSERT_EQ(solution_matrix.get(2, 0), Winterval(43.25, 2845.5));
    ASSERT_EQ(solution_matrix.get(2, 1), Winterval(888.75, 33382));
    ASSERT_EQ(solution_matrix.get(2, 2), Winterval(-2842.5, -40.25));
    ASSERT_EQ(solution_matrix.get(2, 3), Winterval(-33377, -883.75));

    free(initial_conditions);
    initial_conditions = nullptr;
}