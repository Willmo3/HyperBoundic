//
// Created by will on 10/23/25.
//

#include <fstream>

#include "domains/Real.hpp"
#include "gtest/gtest.h"
#include "../../src/solvers/difference/LaxFriedrichsSolver.hpp"
#include "solvers/flux/CubicFlux.hpp"
#include "Waffine/WaffineForm.hpp"
#include "Winterval/Winterval.hpp"
#include "Wixed/WixedForm.hpp"

/**
 * Serialize a real system approximation.
 */
TEST(serialization, serialize_real) {
    // Start with a simple lax friedrichs approximation.

    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_t = 0.02; // Spacing of time. Assuming operating over 4 logical time split into 4 timesteps = 1.
    double delta_x = 1; // Spatial discretization. Assuming total space of four split into 4 parts = 1.

    auto initial_conditions = std::unique_ptr<Real>(static_cast<Real *>(calloc(discretization_size, sizeof(double))));

    initial_conditions.get()[0] = 1.0;
    initial_conditions.get()[1] = 2.0;
    initial_conditions.get()[2] = 3.0;
    initial_conditions.get()[3] = 4.0;

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

    auto initial_conditions = std::unique_ptr<Winterval>(static_cast<Winterval *>(calloc(discretization_size, sizeof(Winterval))));

    initial_conditions.get()[0] = Winterval(0, 1);
    initial_conditions.get()[1] = Winterval(1, 2);
    initial_conditions.get()[2] = Winterval(2, 3);
    initial_conditions.get()[3] = Winterval(3, 4);

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

    auto initial_conditions = std::unique_ptr<WaffineForm>(static_cast<WaffineForm *>(calloc(discretization_size, sizeof(WaffineForm))));

    initial_conditions.get()[0] = WaffineForm(Winterval(0, 1));
    initial_conditions.get()[1] = WaffineForm(Winterval(1, 2));
    initial_conditions.get()[2] = WaffineForm(Winterval(2, 3));
    initial_conditions.get()[3] = WaffineForm(Winterval(3, 4));

    auto solution_matrix = LaxFriedrichsSolver<WaffineForm>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<WaffineForm>());

    // Now, serialize this.

    auto strrep = solution_matrix.to_json_string();

    auto deserialized_matrix = RectangularMesh<WaffineForm>::from_json_string(strrep);
    ASSERT_TRUE(solution_matrix.equals(deserialized_matrix));
}

TEST(serialization, serialize_mixed) {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_t = 0.02;
    double delta_x = 1;

    auto initial_conditions = std::unique_ptr<WixedForm>(static_cast<WixedForm *>(calloc(discretization_size, sizeof(WixedForm))));
    initial_conditions.get()[0] = WixedForm(Winterval(0, 1));
    initial_conditions.get()[1] = WixedForm(Winterval(1, 2));
    initial_conditions.get()[2] = WixedForm(Winterval(2, 3));
    initial_conditions.get()[3] = WixedForm(Winterval(3, 4));

    auto solution_matrix = LaxFriedrichsSolver<WixedForm>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<WixedForm>());

    auto strrep = solution_matrix.to_json_string();

    auto deserialized_matrix = RectangularMesh<WixedForm>::from_json_string(strrep);
    ASSERT_TRUE(solution_matrix.equals(deserialized_matrix));
}