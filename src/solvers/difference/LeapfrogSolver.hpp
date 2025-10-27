//
// Created by will on 10/9/25.
//

#ifndef PDEAPPROX_LEAPFROGSOLVER_H
#define PDEAPPROX_LEAPFROGSOLVER_H
#include <cmath>

#include "LaxFriedrichsSolver.hpp"
#include "domains/Numeric.hpp"
#include "meshes/RectangularMesh.hpp"
#include "../flux/FluxFunction.hpp"
#include "DifferenceHelpers.hpp"

template<typename T>
requires Numeric<T>
class LeapfrogSolver {

public:
    /**
     * @brief Given a set of initial conditions over some discretization of a 1d space, a time discretization, and a number of timesteps,
     * Approximate the values of the system at different points in time using the second-order leapfrog method.
     *
     * @param initial_state Array of size discretization_size representing the initial state of the system
     * @param num_timesteps Number of timesteps for the approximation. (num rows). Must be >= 2 to prime with Lax-Friedrichs.
     * @param discretization_size Size of space being approximated. (num cols)
     * @param delta_t Time discretization... i.e. how much time is passing logically for each step. Must be > 0, < INFINITY
     * @param delta_x Space discretization... i.e. how much space is passing logically for each step. Must be > 0, < INFINITY
     * @param flux Flux function to use for this approximation.
     * @return a discretization of the partial differential equation system.
     */
    static RectangularMesh<T> solve(std::unique_ptr<T> &initial_state, uint32_t discretization_size, uint32_t num_timesteps, double delta_t, double delta_x, FluxFunction<T>* flux) {
        assert(delta_t > 0 && delta_t < INFINITY);
        assert(delta_x > 0 && delta_x < INFINITY);
        assert(num_timesteps >= 2); // Need at least two timesteps to prime with Lax-Friedrichs.

        auto k = delta_t / delta_x;
        auto solution = RectangularMesh<T>(discretization_size, num_timesteps);
        solution.copy_initial_conditions(initial_state);

        auto first_row = LaxFriedrichsSolver<T>::solve(initial_state, discretization_size, 2, delta_t, delta_x, flux);
        // Copy first row of Lax-Friedrichs solution into our solution matrix.
        for (auto x = 0; x < discretization_size; x++) {
            solution.set(1, x, first_row.get(1, x));
        }
        cfl_check_row<T>(solution, flux, delta_t, delta_x, 1);

        // With first timestep primed, move to leapfrog.
        for (auto timestep = 1; timestep < num_timesteps - 1; timestep++) {
            for (auto x = 1; x < discretization_size - 1; x++) {
                auto u_x_plus_1 = solution.get(timestep, x + 1);
                auto u_x_minus_1 = solution.get(timestep, x - 1);
                auto u_x_prev = solution.get(timestep - 1, x);
                solution.set(timestep + 1, x, leapfrog_stencil(u_x_plus_1, u_x_minus_1, u_x_prev, k, flux));
            }
            // Currently, only support periodic boundary conditions.
            solution.set(timestep + 1, 0, leapfrog_stencil(
                solution.get(timestep, 1),
                solution.get(timestep, solution.discretization_size() - 1),
                solution.get(timestep - 1, 0),
                k, flux));

            solution.set(timestep + 1, solution.discretization_size() - 1, leapfrog_stencil(
                solution.get(timestep, 0),
                solution.get(timestep, discretization_size - 2),
                solution.get(timestep - 1, discretization_size - 1),
                k, flux));

            // After each run through, check that CFL satisfied.
            cfl_check_row<T>(solution, flux, delta_t, delta_x, timestep);
        }

        return solution;
    }

    /*
     * Stencils
     */
    static T leapfrog_stencil(T u_x_plus_1, T u_x_minus_1, T u_x_prev, double k, FluxFunction<T> *flux) {
        return u_x_prev - (flux->flux(u_x_plus_1) - flux->flux(u_x_minus_1)) * k;
    }
};

#endif //PDEAPPROX_LEAPFROGSOLVER_H