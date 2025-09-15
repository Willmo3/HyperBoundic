//
// Created by will on 9/11/25.
//

#include <cmath>

#include "LaxFriedrichsSolver.h"

#include <cassert>
#include <functional>

/*
 * Flux functions
 */
double cubic_flux(double value) {
    return pow(value, 3);
}
double unit_flux(double value) {
    return value;
}

/*
 * Stencils
 */
double lax_friedrichs_stencil(double u_i_plus_1, double u_i_minus_1, double k, const std::function<double(double)>& flux) {
    return 0.5 * (u_i_plus_1 + u_i_minus_1) - k * (flux(u_i_plus_1) - flux(u_i_minus_1));
}

/*
 * Approximator
 */
PdeDiscretization lax_friedrichs_solver(const double* initial_state, uint32_t discretization_size, uint32_t num_timesteps, double delta_t, double delta_x) {
    assert(delta_t > 0 && delta_t < INFINITY);
    assert(delta_x > 0 && delta_x < INFINITY);

    auto k = delta_t / delta_x * 1/2;
    auto solution = PdeDiscretization(discretization_size, num_timesteps, initial_state);

    // copy initial state into solution

    for (auto timestep = 0; timestep < num_timesteps - 1; timestep++) {
        for (auto x = 1; x < discretization_size - 1; x++) {
            auto u_x_plus_1 = solution.get(timestep, x + 1);
            auto u_x_minus_1 = solution.get(timestep, x - 1);
            solution.set(timestep + 1, x, lax_friedrichs_stencil(u_x_plus_1, u_x_minus_1, k, unit_flux));
        }

        // Currently, only support periodic boundary conditions.

        solution.set(timestep + 1, 0, lax_friedrichs_stencil(
            solution.get(timestep, 1),
            solution.get(timestep, solution.discretization_size() - 1),
            k, unit_flux));

        solution.set(timestep + 1, solution.discretization_size() - 1, lax_friedrichs_stencil(
            solution.get(timestep, 0),
            solution.get(timestep, discretization_size - 2),
            k, unit_flux));
    }

    return solution;
}