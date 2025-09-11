
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>

#include "src/LaxFriedrichsSolver.h"
#include "src/PdeDiscretization.h"


int main() {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_x = 1; // Spatial discretization. Assuming total space of four split into 4 parts = 1.
    double delta_t = 1; // Spacing of time. Assuming operating over 4 logical time split into 4 timesteps = 1.

    auto initial_conditions = static_cast<double *>(calloc(discretization_size, sizeof(double)));
    assert(initial_conditions);

    initial_conditions[0] = 1.0;
    initial_conditions[1] = 2.0;
    initial_conditions[2] = 3.0;
    initial_conditions[3] = 4.0;

    auto solution_matrix = lax_friedrichs_solver(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x);
    solution_matrix.print_system();

    free(initial_conditions);
    initial_conditions = nullptr;
    return 0;
}