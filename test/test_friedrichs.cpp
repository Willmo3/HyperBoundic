//
// Created by will on 9/18/25.
//

#include <gtest/gtest.h>
#include "../src/domains/Real.hpp"
#include "../src/solvers/LaxFriedrichsSolver.hpp"
#include "../src/solvers/flux/CubicFlux.hpp"
#include "../lib/Winterval/src/Winterval.h"

void assert_eq_bounded_interval(Winterval a, Winterval b) {
    ASSERT_NEAR(a.min(), b.min(), 1e-5);
    ASSERT_NEAR(a.max(), b.max(), 1e-5);
}

// Note: with cubic flux, delta_t must be very low -- this function is highly vulnerable to blowup!
TEST(friedrichs, real_approx) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_t = 0.02; // Spacing of time. Assuming operating over 4 logical time split into 4 timesteps = 1.
    double delta_x = 1; // Spatial discretization. Assuming total space of four split into 4 parts = 1.

    auto initial_conditions = static_cast<Real *>(calloc(discretization_size, sizeof(double)));
    assert(initial_conditions);

    initial_conditions[0] = 1.0;
    initial_conditions[1] = 2.0;
    initial_conditions[2] = 3.0;
    initial_conditions[3] = 4.0;

    auto solution_matrix = LaxFriedrichsSolver<Real>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<Real>());
    // test the third row of the discretization.
    ASSERT_NEAR(solution_matrix.get(2, 0).value(), 2.062752, 0.001);
    ASSERT_NEAR(solution_matrix.get(2, 1).value(), 3.305912, 0.001);
    ASSERT_NEAR(solution_matrix.get(2, 2).value(), 1.937248, 0.001);
    ASSERT_NEAR(solution_matrix.get(2, 3).value(), 2.694088, 0.001);

    free(initial_conditions);
    initial_conditions = nullptr;
}

TEST(friedrichs, interval_approx) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_t = 0.02;
    double delta_x = 1;

    auto initial_conditions = static_cast<Winterval *>(calloc(discretization_size, sizeof(Winterval)));
    assert(initial_conditions);

    initial_conditions[0] = Winterval(0, 1);
    initial_conditions[1] = Winterval(1, 2);
    initial_conditions[2] = Winterval(2, 3);
    initial_conditions[3] = Winterval(3, 4);

    auto solution_matrix = LaxFriedrichsSolver<Winterval>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<Winterval>());
    assert_eq_bounded_interval(solution_matrix.get(2, 0), Winterval(0.840360, 2.213081));
    assert_eq_bounded_interval(solution_matrix.get(2, 1), Winterval(1.663154, 3.672608));
    assert_eq_bounded_interval(solution_matrix.get(2, 2), Winterval(0.786919, 2.159640));
    assert_eq_bounded_interval(solution_matrix.get(2, 3),  Winterval(1.327392, 3.336846));

    free(initial_conditions);
    initial_conditions = nullptr;
}