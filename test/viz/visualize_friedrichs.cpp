//
// Created by will on 10/17/25.
//

#include "solvers/difference/LaxFriedrichsSolver.hpp"
#include "domains/Real.hpp"
#include "flux/CubicFlux.hpp"
#include "visualization/MeshVisualizer.hpp"

int main() {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_x = 1; // Spatial discretization. Assuming total space of four split into 4 parts = 1.
    double delta_t = 0.02; // Spacing of time. Assuming operating over 4 logical time split into 4 timesteps = 1.

    auto initial_conditions = std::vector<Real>(discretization_size);

    initial_conditions[0] = 1.0;
    initial_conditions[1] = 2.0;
    initial_conditions[2] = 3.0;
    initial_conditions[3] = 4.0;

    auto solution_matrix = LaxFriedrichsSolver(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<Real>()).solve();
    solution_matrix.print_system();
    show_real_surface(&solution_matrix);

    return 0;
}