//
// Created by will on 10/31/25.
//

// Sanity tests for various flux functions.

#include <gtest/gtest.h>

#include <cstdint>
#include <vector>

#include "domains/Real.hpp"
#include "solvers/difference/LaxFriedrichsSolver.hpp"
#include "flux/BuckleyLeverettFlux.hpp"
#include "flux/BurgersFlux.hpp"
#include "flux/LwrFlux.hpp"

RectangularMesh<Real> solve_flux(FluxFunction<Real> * f) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 5;
    double delta_t = 0.02; // Spacing of time. Assuming operating over 4 logical time split into 4 timesteps = 1.
    double delta_x = 1; // Spatial discretization. Assuming total space of four split into 4 parts = 1.

    auto initial_conditions = std::vector<Real>(discretization_size);

    initial_conditions[0] = 1.0;
    initial_conditions[1] = 2.0;
    initial_conditions[2] = 3.0;
    initial_conditions[3] = 4.0;

    return LaxFriedrichsSolver<Real>().solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, f);
}

TEST(flux, burgers_flux) {
    auto solution_matrix = solve_flux(new BurgersFlux<Real>());
    ASSERT_NEAR(solution_matrix.get(4, 0).value(), 1.999997, 0.000001);
    ASSERT_NEAR(solution_matrix.get(4, 1).value(), 2.999987, 0.000001);
    ASSERT_NEAR(solution_matrix.get(4, 2).value(), 2.000003, 0.000001);
    ASSERT_NEAR(solution_matrix.get(4, 3).value(), 3.000013, 0.000001);
}

TEST(flux, lwr_flux) {
    auto solution_matrix = solve_flux(new LwrFlux<Real>());
    ASSERT_NEAR(solution_matrix.get(4, 0).value(), 1.999987, 0.000001);
    ASSERT_NEAR(solution_matrix.get(4, 1).value(), 2.999900, 0.000001);
    ASSERT_NEAR(solution_matrix.get(4, 2).value(), 2.000013, 0.000001);
    ASSERT_NEAR(solution_matrix.get(4, 3).value(), 3.000100, 0.000001);
}

TEST(flux, buckley_lev_flux) {
    auto solution_matrix = solve_flux(new BuckleyLeverett<Real>());
    ASSERT_NEAR(solution_matrix.get(2, 0).value(), 2.000001, 0.000001);
    ASSERT_NEAR(solution_matrix.get(2, 1).value(), 3.000000, 0.000001);
    ASSERT_NEAR(solution_matrix.get(2, 2).value(), 1.999999, 0.000001);
    ASSERT_NEAR(solution_matrix.get(2, 3).value(), 3.000000, 0.000001);
}