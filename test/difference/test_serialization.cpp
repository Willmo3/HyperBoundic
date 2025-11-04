//
// Created by will on 10/23/25.
//

#include <fstream>

#include "domains/Real.hpp"
#include "gtest/gtest.h"
#include "solvers/difference/LaxFriedrichsSolver.hpp"
#include "flux/CubicFlux.hpp"
#include "Caffeine/AffineForm.hpp"
#include "Winterval/Winterval.hpp"
#include "DualDomain/MixedForm.hpp"

/**
 * Serialize a real system approximation.
 */
TEST(serialization, serialize_real) {
    // Start with a simple lax friedrichs approximation.

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

    // Now, serialize this.

    auto strrep = solution_matrix.to_json_string();

    auto deserialized_matrix = RectangularMesh<Real>::from_json_string(strrep);
    ASSERT_TRUE(solution_matrix.equals(deserialized_matrix));
}

TEST(serialization, serialize_interval) {
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

    // Now, serialize this.

    auto strrep = solution_matrix.to_json_string();

    auto deserialized_matrix = RectangularMesh<Winterval>::from_json_string(strrep);
    ASSERT_TRUE(solution_matrix.equals(deserialized_matrix));
}

TEST(serialization, serialize_affine) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_x = 1;
    double delta_t = 0.02;

    auto initial_conditions = std::vector<AffineForm>(discretization_size);

    initial_conditions[0] = AffineForm(Winterval(0, 1));
    initial_conditions[1] = AffineForm(Winterval(1, 2));
    initial_conditions[2] = AffineForm(Winterval(2, 3));
    initial_conditions[3] = AffineForm(Winterval(3, 4));

    auto solution_matrix = LaxFriedrichsSolver<AffineForm>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<AffineForm>());

    // Now, serialize this.

    auto strrep = solution_matrix.to_json_string();

    auto deserialized_matrix = RectangularMesh<AffineForm>::from_json_string(strrep);
    ASSERT_TRUE(solution_matrix.equals(deserialized_matrix));
}

TEST(serialization, serialize_mixed) {
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

    auto strrep = solution_matrix.to_json_string();

    auto deserialized_matrix = RectangularMesh<MixedForm>::from_json_string(strrep);
    ASSERT_TRUE(solution_matrix.equals(deserialized_matrix));
}