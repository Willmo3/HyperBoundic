
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>

/*
 * Stencils
 */
double lax_friedrichs_stencil_cubic(double u_i_plus_1, double u_i_minus_1, double k) {
    return 0.5 * (u_i_plus_1 + u_i_minus_1) - k * (pow(u_i_plus_1, 3) - pow(u_i_minus_1, 3));
}
double lax_friedrichs_stencil_simple(double u_i_plus_1, double u_i_minus_1, double k) {
    return 0.5 * (u_i_plus_1 + u_i_minus_1) - k * (u_i_plus_1 - u_i_minus_1);
}

void print_system(double *system, uint32_t discretization_size, uint32_t num_timesteps) {
    for (auto t = 0; t < num_timesteps; t++) {
        std::cout << "T" << t << ": ";
        for (auto i = 0; i < discretization_size; i++) {
            std::cout << system[discretization_size * t + i] << " ";
        }
        std::cout << std::endl;
    }
}

// Insight: conservation law formulation, hence = 0 -- very convenient!
//  -- sum energy in system cannot change.
double *lax_friedrichs_solver(double *initial_state, uint32_t discretization_size, uint32_t num_timesteps, double delta_t, double delta_x) {
    auto k = delta_t / delta_x * 1/2;
    auto solution = static_cast<double *>(calloc(sizeof(double), num_timesteps * discretization_size));
    assert(solution);

    // Copy initial values into solution matrix
    assert(memcpy(solution, initial_state, discretization_size * sizeof(double)));

    // copy initial state into solution

    for (auto timestep = 0; timestep < num_timesteps - 1; timestep++) {
        for (auto x = 1; x < discretization_size - 1; x++) {
            auto u_x_plus_1 = solution[timestep * discretization_size + x + 1];
            auto u_x_minus_1 = solution[timestep * discretization_size + x - 1];

            solution[(timestep + 1) * discretization_size + x] =
                lax_friedrichs_stencil_simple(u_x_plus_1, u_x_minus_1, k);
        }

        // Currently, only support periodic boundary conditions.
        solution[(timestep + 1) * discretization_size] = lax_friedrichs_stencil_simple(
            solution[timestep * discretization_size + 1],
            solution[timestep * discretization_size + discretization_size - 1],
            k);

        solution[(timestep + 1) * discretization_size + discretization_size - 1] = lax_friedrichs_stencil_simple(
            solution[timestep * discretization_size],
            solution[timestep * discretization_size + discretization_size - 2],
            k);
    }

    return solution;
}

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
    print_system(solution_matrix, discretization_size, num_timesteps);

    free(initial_conditions);
    initial_conditions = nullptr;
    free(solution_matrix);
    solution_matrix = nullptr;
    assert(!initial_conditions);
    assert(!solution_matrix);
    return 0;
}


// alternative: define k elsewhere.

// Insight: change in dependent variable wrt time, plus change in flux spatialy, equals zero.
// Solving this PDE shows us the dependent variable at a given space and time.

// Hyperbolic PDE: real, distinct eigenvalues representing propagation speed.

