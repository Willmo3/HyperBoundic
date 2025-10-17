//
// Created by will on 10/17/25.
//


#include <cassert>
#include <cstdint>
#include <memory>
#include "../../src/solvers/LeapfrogSolver.hpp"
#include "../../src/domains/Real.hpp"
#include "../../src/solvers/flux/CubicFlux.hpp"
#include "../../src/visualization/DiscretizationVisualizers.hpp"

int main() {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_x = 1; // Spatial discretization. Assuming total space of four split into 4 parts = 1.
    double delta_t = 0.02; // Spacing of time. Assuming operating over 4 logical time split into 4 timesteps = 1.

    auto initial_conditions = std::unique_ptr<Real>(static_cast<Real *>(calloc(discretization_size, sizeof(double))));
    assert(initial_conditions);

    initial_conditions.get()[0] = 1.0;
    initial_conditions.get()[1] = 2.0;
    initial_conditions.get()[2] = 3.0;
    initial_conditions.get()[3] = 4.0;

    auto solution_matrix = LeapfrogSolver<Real>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<Real>());
    solution_matrix.print_system();
    show_real_surface(&solution_matrix);
}
