//
// Created by will on 9/11/25.
//

#include <cmath>

#include "LaxFriedrichsSolver.h"
/*
 * Stencils
 */
double lax_friedrichs_stencil_cubic(double u_i_plus_1, double u_i_minus_1, double k) {
    return 0.5 * (u_i_plus_1 + u_i_minus_1) - k * (pow(u_i_plus_1, 3) - pow(u_i_minus_1, 3));
}
double lax_friedrichs_stencil_simple(double u_i_plus_1, double u_i_minus_1, double k) {
    return 0.5 * (u_i_plus_1 + u_i_minus_1) - k * (u_i_plus_1 - u_i_minus_1);
}

// Insight: conservation law formulation, hence = 0 -- very convenient!
//  -- sum energy in system cannot change.
PdeDiscretization lax_friedrichs_solver(const double* initial_state, uint32_t discretization_size, uint32_t num_timesteps, double delta_t, double delta_x) {
    auto k = delta_t / delta_x * 1/2;
    auto solution = PdeDiscretization(discretization_size, num_timesteps, initial_state);

    // copy initial state into solution

    for (auto timestep = 0; timestep < num_timesteps - 1; timestep++) {
        for (auto x = 1; x < discretization_size - 1; x++) {
            auto u_x_plus_1 = solution.get(timestep, x + 1);
            auto u_x_minus_1 = solution.get(timestep, x - 1);
            solution.set(timestep + 1, x, lax_friedrichs_stencil_simple(u_x_plus_1, u_x_minus_1, k));
        }

        // Currently, only support periodic boundary conditions.

        solution.set(timestep + 1, 0, lax_friedrichs_stencil_simple(
            solution.get(timestep, 1),
            solution.get(timestep, solution.discretization_size() - 1),
            k));

        solution.set(timestep + 1, solution.discretization_size() - 1, lax_friedrichs_stencil_simple(
            solution.get(timestep, 0),
            solution.get(timestep, discretization_size - 2),
            k));
    }

    return solution;
}