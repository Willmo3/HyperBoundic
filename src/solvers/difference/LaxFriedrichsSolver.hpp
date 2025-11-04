//
// Created by will on 9/11/25.
//

#ifndef PDENCLOSE_LAXFRIEDRICHSSOLVER_H
#define PDENCLOSE_LAXFRIEDRICHSSOLVER_H

#include <cmath>

#include "domains/Numeric.hpp"
#include "meshes/RectangularMesh.hpp"
#include "flux/FluxFunction.hpp"
#include "DifferenceHelpers.hpp"

// Using template classes to ease instantiation.
template<typename T>
requires Numeric<T>
class LaxFriedrichsSolver {

public:
    /**
     * @brief Given a set of initial conditions over some discretization of a 1d space, a time discretization, and a number of timesteps,
     * Approximate the values of the system at different points in time.
     *
     * @param initial_state vector of size discretization_size representing the initial state of the system
     * @param num_timesteps Number of timesteps for the approximation. (num rows)
     * @param discretization_size Size of space being approximated. (num cols)
     * @param delta_t Time discretization... i.e. how much time is passing logically for each step. Must be > 0, < INFINITY
     * @param delta_x Space discretization... i.e. how much space is passing logically for each step. Must be > 0, < INFINITY
     * @param flux Flux function to use for this approximation.
     * @return a discretization of the partial differential equation system.
     */
    static RectangularMesh<T> solve(const std::vector<T> &initial_state, uint32_t discretization_size, uint32_t num_timesteps, double delta_t, double delta_x, FluxFunction<T>* flux) {
        assert(delta_t > 0 && delta_t < INFINITY);
        assert(delta_x > 0 && delta_x < INFINITY);

        auto k = delta_t / delta_x * 1/2;
        auto solution = RectangularMesh<T>(discretization_size, num_timesteps);
        solution.copy_initial_conditions(initial_state);

        for (auto timestep = 0; timestep < num_timesteps - 1; timestep++) {
            for (auto x = 1; x < discretization_size - 1; x++) {
                auto u_x_plus_1 = solution.get(timestep, x + 1);
                auto u_x_minus_1 = solution.get(timestep, x - 1);
                solution.set(timestep + 1, x, lax_friedrichs_stencil(u_x_plus_1, u_x_minus_1, k, flux));
            }

            // Currently, only support periodic boundary conditions.
            solution.set(timestep + 1, 0, lax_friedrichs_stencil(
                solution.get(timestep, 1),
                solution.get(timestep, solution.discretization_size() - 1),
                k, flux));

            solution.set(timestep + 1, solution.discretization_size() - 1, lax_friedrichs_stencil(
                solution.get(timestep, 0),
                solution.get(timestep, discretization_size - 2),
                k, flux));

            // After each run through, check that CFL satisfied.
            cfl_check_row<T>(solution, flux, delta_t, delta_x, timestep);
        }

        return solution;
    }

    /*
     * Stencils
     */
    static T lax_friedrichs_stencil(T u_i_plus_1, T u_i_minus_1, double k, FluxFunction<T> *flux) {
        return (u_i_plus_1 + u_i_minus_1) * 0.5 - (flux->flux(u_i_plus_1) - flux->flux(u_i_minus_1)) * k;
    }
};

#endif //PDENCLOSE_LAXFRIEDRICHSSOLVER_H