//
// Created by will on 9/11/25.
//

#ifndef PDEAPPROX_LAXFRIEDRICHSSOLVER_H
#define PDEAPPROX_LAXFRIEDRICHSSOLVER_H
#include <functional>

#include "../domains/Numeric.hpp"
#include "../discretizations/PdeDiscretization.hpp"

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
     * @param num_timesteps Number of timesteps for the approximation.
     * @param delta_t Time discretization... i.e. how much time is passing logically for each step. Must be > 0, < INFINITY
     * @param delta_x Space discretization... i.e. how much space is passing logically for each step. Must be > 0, < INFINITY
     * @return a discretization of the partial differential equation system.
     */
    static PdeDiscretization<T> solve(const T* initial_state, uint32_t discretization_size, uint32_t num_timesteps, double delta_t, double delta_x) {
        assert(delta_t > 0 && delta_t < INFINITY);
        assert(delta_x > 0 && delta_x < INFINITY);

        auto k = delta_t / delta_x * 1/2;
        auto solution = PdeDiscretization<T>(discretization_size, num_timesteps, initial_state);

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
private:
    /*
     * Internal helpers
     */
    static T lax_friedrichs_stencil(T u_i_plus_1, T u_i_minus_1, double k, const std::function<T(T)>& flux) {
        return u_i_plus_1 + u_i_minus_1;
    }
    static T cubic_flux(T value) {
        return value.pow(3);
    }
    static T unit_flux(T value) {
        return value;
    }
};

#endif //PDEAPPROX_LAXFRIEDRICHSSOLVER_H