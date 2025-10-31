//
// Created by will on 9/18/25.
//

#include <gtest/gtest.h>

#include "Caffeine/AffineForm.hpp"
#include "domains/Real.hpp"
#include "../../src/solvers/difference/LaxFriedrichsSolver.hpp"
#include "solvers/flux/CubicFlux.hpp"
#include "Winterval/Winterval.hpp"
#include "DualDomain/MixedForm.hpp"
#include "solvers/flux/BurgersFlux.hpp"

void assert_eq_bounded_interval(Winterval a, Winterval b) {
    ASSERT_NEAR(a.min(), b.min(), 1e-5);
    ASSERT_NEAR(a.max(), b.max(), 1e-5);
}

// TODO: tests weak, not wide enough -- too many boundary conditions

// Note: with cubic flux, delta_t must be very low -- this function is highly vulnerable to blowup!
TEST(friedrichs, real_approx_cubic_flux) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_t = 0.02; // Spacing of time. Assuming operating over 4 logical time split into 4 timesteps = 1.
    double delta_x = 1; // Spatial discretization. Assuming total space of four split into 4 parts = 1.

    auto initial_conditions = std::vector<Real>(discretization_size);

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
}

TEST(friedrichs, real_approx_burgers_flux) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 5;
    double delta_t = 0.02; // Spacing of time. Assuming operating over 4 logical time split into 4 timesteps = 1.
    double delta_x = 1; // Spatial discretization. Assuming total space of four split into 4 parts = 1.

    auto initial_conditions = std::vector<Real>(discretization_size);

    initial_conditions[0] = 1.0;
    initial_conditions[1] = 2.0;
    initial_conditions[2] = 3.0;
    initial_conditions[3] = 4.0;

    auto solution_matrix = LaxFriedrichsSolver<Real>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new BurgersFlux<Real>());
    ASSERT_NEAR(solution_matrix.get(4, 0).value(), 1.999997, 0.000001);
    ASSERT_NEAR(solution_matrix.get(4, 1).value(), 2.999987, 0.000001);
    ASSERT_NEAR(solution_matrix.get(4, 2).value(), 2.000003, 0.000001);
    ASSERT_NEAR(solution_matrix.get(4, 3).value(), 3.000013, 0.000001);
}

TEST(friedrichs, interval_approx) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_t = 0.02;
    double delta_x = 1;

    auto initial_conditions = std::vector<Winterval>(discretization_size);

    initial_conditions[0] = Winterval(0, 1);
    initial_conditions[1] = Winterval(1, 2);
    initial_conditions[2] = Winterval(2, 3);
    initial_conditions[3] = Winterval(3, 4);

    auto solution_matrix = LaxFriedrichsSolver<Winterval>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<Winterval>());
    assert_eq_bounded_interval(solution_matrix.get(2, 0), Winterval(0.840360, 2.213081));
    assert_eq_bounded_interval(solution_matrix.get(2, 1), Winterval(1.663154, 3.672608));
    assert_eq_bounded_interval(solution_matrix.get(2, 2), Winterval(0.786919, 2.159640));
    assert_eq_bounded_interval(solution_matrix.get(2, 3),  Winterval(1.327392, 3.336846));
}

TEST(friedrichs, affine_approx) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_t = 0.02;
    double delta_x = 1;

    auto initial_conditions = std::vector<AffineForm>(discretization_size);

    initial_conditions[0] = AffineForm(Winterval(0, 1));
    initial_conditions[1] = AffineForm(Winterval(1, 2));
    initial_conditions[2] = AffineForm(Winterval(2, 3));
    initial_conditions[3] = AffineForm(Winterval(3, 4));

    auto solution_matrix = LaxFriedrichsSolver<AffineForm>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<AffineForm>());
    assert_eq_bounded_interval(solution_matrix.get(2, 0).to_interval(), Winterval(0.939509, 2.102490));
    assert_eq_bounded_interval(solution_matrix.get(2, 1).to_interval(), Winterval(1.932881, 3.365835));
    assert_eq_bounded_interval(solution_matrix.get(2, 2).to_interval(), Winterval(0.951364, 2.006637));
    assert_eq_bounded_interval(solution_matrix.get(2, 3).to_interval(), Winterval(1.877454, 2.823831));
}

TEST(friedrichs, mixed_approx) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_t = 0.02;
    double delta_x = 1;

    auto initial_conditions = std::vector<MixedForm>(discretization_size);
    initial_conditions[0] = MixedForm(Winterval(0, 1));
    initial_conditions[1] = MixedForm(Winterval(1, 2));
    initial_conditions[2] = MixedForm(Winterval(2, 3));
    initial_conditions[3] = MixedForm(Winterval(3, 4));

    auto solution_matrix = LaxFriedrichsSolver<MixedForm>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<MixedForm>());
    // Notice that mixed approximation is slighly tighter than pure affine at index (2,0).
    assert_eq_bounded_interval(solution_matrix.get(2, 0).interval_bounds(), Winterval(0.951795, 2.102490));
    assert_eq_bounded_interval(solution_matrix.get(2, 1).interval_bounds(), Winterval(1.932881, 3.365835));
    assert_eq_bounded_interval(solution_matrix.get(2, 2).interval_bounds(), Winterval(0.951364, 2.006637));
    assert_eq_bounded_interval(solution_matrix.get(2, 3).interval_bounds(), Winterval(1.877454, 2.823831));
}