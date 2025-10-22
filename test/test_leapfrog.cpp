//
// Created by will on 10/17/25.
//

#include "Waffine/WaffineForm.hpp"
#include "../src/domains/Real.hpp"
#include "../src/solvers/LeapfrogSolver.hpp"
#include "../src/solvers/flux/CubicFlux.hpp"
#include "Winterval/Winterval.hpp"
#include "Wixed/WixedForm.hpp"
#include "gtest/gtest.h"

void assert_eq_bounded_interval(Winterval a, Winterval b) {
    ASSERT_NEAR(a.min(), b.min(), 1e-5);
    ASSERT_NEAR(a.max(), b.max(), 1e-5);
}

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

TEST(leapfrog, interval_approx) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_t = 0.02;
    double delta_x = 1;

    auto initial_conditions = std::unique_ptr<Winterval>(static_cast<Winterval *>(calloc(discretization_size, sizeof(Winterval))));

    initial_conditions.get()[0] = Winterval(0, 1);
    initial_conditions.get()[1] = Winterval(1, 2);
    initial_conditions.get()[2] = Winterval(2, 3);
    initial_conditions.get()[3] = Winterval(3, 4);

    auto solution_matrix = LeapfrogSolver<Winterval>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<Winterval>());
    assert_eq_bounded_interval(Winterval(-0.119280, 1.226161), solution_matrix.get(2, 0));
    assert_eq_bounded_interval(Winterval(0.766308, 2.905216), solution_matrix.get(2, 1));
    assert_eq_bounded_interval(Winterval(1.773839, 3.119280), solution_matrix.get(2, 2));
    assert_eq_bounded_interval(Winterval(2.094784, 4.233692), solution_matrix.get(2, 3));
}

TEST(leapfrog, affine_approx) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_x = 1;
    double delta_t = 0.02;

    auto initial_conditions = std::unique_ptr<WaffineForm>(static_cast<WaffineForm *>(calloc(discretization_size, sizeof(WaffineForm))));

    initial_conditions.get()[0] = WaffineForm(Winterval(0, 1));
    initial_conditions.get()[1] = WaffineForm(Winterval(1, 2));
    initial_conditions.get()[2] = WaffineForm(Winterval(2, 3));
    initial_conditions.get()[3] = WaffineForm(Winterval(3, 4));

    auto solution_matrix = LeapfrogSolver<WaffineForm>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<WaffineForm>());
    assert_eq_bounded_interval(solution_matrix.get(2, 0).to_interval(), Winterval(-0.076409, 1.160407));
    assert_eq_bounded_interval(solution_matrix.get(2, 1).to_interval(), Winterval(0.924492, 2.672938));
    assert_eq_bounded_interval(solution_matrix.get(2, 2).to_interval(), Winterval(1.918658, 2.997344));
    assert_eq_bounded_interval(solution_matrix.get(2, 3).to_interval(), Winterval(2.728068, 3.674502));
}

TEST(leapfrog, affine_approx_mixed) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_x = 1;
    double delta_t = 0.02;

    auto initial_conditions = std::unique_ptr<WixedForm>(static_cast<WixedForm *>(calloc(discretization_size, sizeof(WixedForm))));

    initial_conditions.get()[0] = WixedForm(Winterval(0, 1));
    initial_conditions.get()[1] = WixedForm(Winterval(1, 2));
    initial_conditions.get()[2] = WixedForm(Winterval(2, 3));
    initial_conditions.get()[3] = WixedForm(Winterval(3, 4));

    auto solution_matrix = LeapfrogSolver<WixedForm>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<WixedForm>());
    assert_eq_bounded_interval(solution_matrix.get(2, 0).interval_bounds(), Winterval(-0.076409,1.160407));
    assert_eq_bounded_interval(solution_matrix.get(2, 1).interval_bounds(), Winterval(0.924492, 2.672938));
    assert_eq_bounded_interval(solution_matrix.get(2, 2).interval_bounds(), Winterval(1.918658, 2.997344));
    assert_eq_bounded_interval(solution_matrix.get(2, 3).interval_bounds(), Winterval(2.728068, 3.674502));
}