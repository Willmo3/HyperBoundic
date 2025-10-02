//
// Created by will on 9/11/25.
//

#ifndef PDEAPPROX_LAXFRIEDRICHSSOLVER_H
#define PDEAPPROX_LAXFRIEDRICHSSOLVER_H

#include <cmath>

#include "../domains/Numeric.hpp"
#include "../discretizations/PdeDiscretization.hpp"
#include "cfl/CflCheck.hpp"

#include "flux/FluxFunction.hpp"

// Using template classes to ease instantiation.
template<typename T>
requires Numeric<T>
class LaxFriedrichsSolver {

public:
    /**
     * @brief Given a set of initial conditions over some discretization of a 1d space, a time discretization, and a number of timesteps,
     * Approximate the values of the system at different points in time.
     *
     * @param initial_state Array of size discretization_size representing the initial state of the system
     * @param num_timesteps Number of timesteps for the approximation. (num rows)
     * @param discretization_size Size of space being approximated. (num cols)
     * @param delta_t Time discretization... i.e. how much time is passing logically for each step. Must be > 0, < INFINITY
     * @param delta_x Space discretization... i.e. how much space is passing logically for each step. Must be > 0, < INFINITY
     * @param flux Flux function to use for this approximation.
     * @return a discretization of the partial differential equation system.
     */
    static PdeDiscretization<T> solve(std::unique_ptr<T> &initial_state, uint32_t discretization_size, uint32_t num_timesteps, double delta_t, double delta_x, FluxFunction<T>* flux) {
        assert(delta_t > 0 && delta_t < INFINITY);
        assert(delta_x > 0 && delta_x < INFINITY);

        auto k = delta_t / delta_x * 1/2;
        auto solution = PdeDiscretization<T>(discretization_size, num_timesteps, initial_state);

        // CFL condition check
        // at each timestep, ensure condition met.
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
            for (auto point = 0; point < discretization_size; point++) {
                if (!cfl_check(flux, solution.get(timestep, point), delta_t, delta_x)) {
                    std::cerr << "System blowup detected, exiting." << std::endl;
                    exit(1);
                }
            }
        }

        return solution;
    }
private:
    /*
     * Internal helpers
     */
    static T lax_friedrichs_stencil(T u_i_plus_1, T u_i_minus_1, double k, FluxFunction<T> *flux) {
        // Keeping these as separate variables allows more complex domains to be operated on as references, rather than values.
        return (u_i_plus_1 + u_i_minus_1) * 0.5 - (flux->flux(u_i_plus_1) - flux->flux(u_i_minus_1)) * k;
    }
};

#endif //PDEAPPROX_LAXFRIEDRICHSSOLVER_H