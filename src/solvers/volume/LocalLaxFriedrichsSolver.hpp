//
// Created by will on 10/27/25.
//

#ifndef PDEAPPROX_LOCALLAXFRIEDRICHSSOLVER_H
#define PDEAPPROX_LOCALLAXFRIEDRICHSSOLVER_H
#include <cmath>

#include "meshes/RectangularMesh.hpp"

template<typename T>
requires Numeric<T>
class LocalLaxFriedrichsSolver {
public:
    /**
     * @brief Approximate a finite volume mesh of a discretized system using the local lax friedrichs solver.
     * The key difference from the standard LF method is that the mesh may be irregular, so the discretization constant
     * is calculated for each control volume cell.
     *
     * @param initial_state Initial values of system.
     * @param width_values Width for each point in the initial state of the system -- irregular mesh.
     * This subsumes delta_x in the standard lax friedrichs solver.
     * @param discretization_size Number of control volume cells.
     * @param num_timesteps Number of timesteps for simulation.
     * @param delta_t Change in time at each point.
     * @param flux Flux function to approximate system with.
     * @return The approximation of the system.
     */
    static RectangularMesh<T> solve(const std::vector<T> &initial_state, const std::vector<double> &width_values, uint32_t discretization_size,
                                    uint32_t num_timesteps, double delta_t, FluxFunction<T> *flux) {
        assert(delta_t > 0 && delta_t < INFINITY);
        // Each point in the discretization must have a corresponding delta_x
        assert(width_values.size() == discretization_size);
        for (auto width_value : width_values) {
            // Verify delta_x for each width.
            assert(width_value > 0 && width_value < INFINITY);
        }

        auto solution = RectangularMesh<T>(discretization_size, num_timesteps);
        solution.copy_initial_conditions(initial_state);

        for (auto timestep = 0; timestep < num_timesteps - 1; timestep++) {
            for (auto x = 1; x < discretization_size - 1; x++) {
                auto k = delta_t / width_values[x] * 1/2;

                auto u_x_plus_1 = solution.get(timestep, x + 1);
                auto u_x_minus_1 = solution.get(timestep, x - 1);
                solution.set(timestep + 1, x, local_lax_friedrichs_stencil(u_x_plus_1, u_x_minus_1, k, flux));
            }

            auto k_0 = delta_t / width_values[0] * 1/2;
            // Currently, only support periodic boundary conditions.
            solution.set(timestep + 1, 0, local_lax_friedrichs_stencil(
                solution.get(timestep, 1),
                solution.get(timestep, solution.discretization_size() - 1),
                k_0, flux));

            auto k_n = delta_t / width_values[solution.discretization_size() - 1] * 1/2;
            solution.set(timestep + 1, solution.discretization_size() - 1, local_lax_friedrichs_stencil(
                solution.get(timestep, 0),
                solution.get(timestep, discretization_size - 2),
                k_n, flux));

            // After each run through, check that CFL satisfied.
            // TODO: cfl_check_row_irregular?
            for (auto x = 0; x < solution.discretization_size(); x++) {
                cfl_check(flux, solution.get(timestep, x), delta_t, width_values[x]);
            }
        }

        return solution;
    }

    // LLF stenci stencil derived from: https://diogenes.bg/ijam/contents/2022-35-4/7/7.pdf
    // Control flow is the same as Lax Friedrichs. Big Rusanov contribution is the flux fn, which we haven't yet implemented.
    // This stencil imposes the discretization on the flux.
    static T local_lax_friedrichs_stencil(T u_i_plus_1, T u_i_minus_1, double k, FluxFunction<T> *flux) {
        return (u_i_plus_1 + u_i_minus_1) * 0.5 - (flux->flux(u_i_plus_1) - flux->flux(u_i_minus_1)) * k * 0.5;
    }
};

#endif //PDEAPPROX_LOCALLAXFRIEDRICHSSOLVER_H