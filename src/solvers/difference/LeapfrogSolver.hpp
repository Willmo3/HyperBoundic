//
// Created by will on 10/9/25.
//

#ifndef PDENCLOSE_LEAPFROGSOLVER_H
#define PDENCLOSE_LEAPFROGSOLVER_H
#include <cmath>

#include "LaxFriedrichsSolver.hpp"
#include "domains/Numeric.hpp"
#include "meshes/Rectangularmesh.hpp"
#include "flux/FluxFunction.hpp"

template<typename T>
requires Numeric<T>
class LeapfrogSolver final: public DifferenceSolver<T> {

public:
    /*
     * Constructor
     */
    LeapfrogSolver(const std::vector<T> &initial_state,
                        uint32_t discretization_size,
                        uint32_t num_timesteps,
                        double delta_t,
                        double delta_x,
                        FluxFunction<T> *flux):

        DifferenceSolver<T>(initial_state, discretization_size, num_timesteps, delta_t, delta_x, flux) {
            assert(num_timesteps >= 2); // Need at least two timesteps to prime with Lax-Friedrichs.
        }

    /*
     * Super fields
     */
    using DifferenceSolver<T>::delta_t;
    using DifferenceSolver<T>::delta_x;
    using DifferenceSolver<T>::mesh;
    using DifferenceSolver<T>::flux;

    void solve() override {
        auto k = delta_t / delta_x;

        // bring back initial state from mesh.
        auto initial_state = std::vector<T>(mesh.discretization_size());
        for (auto x = 0; x < mesh.discretization_size(); x++) {
            initial_state[x] = mesh.get(0, x);
        }

        auto first_row = LaxFriedrichsSolver<T>(initial_state, mesh.discretization_size(), 2, delta_t, delta_x, flux).solve();
        // Copy first row of Lax-Friedrichs solution into our solution matrix.
        for (auto x = 0; x < mesh.discretization_size(); x++) {
            mesh.set(1, x, first_row.get(1, x));
        }
       this->cfl_check_row(mesh, flux, delta_t, delta_x, 1);

        // With first timestep primed, move to leapfrog.
        for (auto timestep = 1; timestep < mesh.num_timesteps() - 1; timestep++) {
            for (auto x = 1; x < mesh.discretization_size() - 1; x++) {
                auto u_x_plus_1 = mesh.get(timestep, x + 1);
                auto u_x_minus_1 = mesh.get(timestep, x - 1);
                auto u_x_prev = mesh.get(timestep - 1, x);
                mesh.set(timestep + 1, x, leapfrog_stencil(u_x_plus_1, u_x_minus_1, u_x_prev, k, flux));
            }
            // Currently, only support periodic boundary conditions.
            mesh.set(timestep + 1, 0, leapfrog_stencil(
                mesh.get(timestep, 1),
                mesh.get(timestep, mesh.discretization_size() - 1),
                mesh.get(timestep - 1, 0),
                k, flux));

            mesh.set(timestep + 1, mesh.discretization_size() - 1, leapfrog_stencil(
                mesh.get(timestep, 0),
                mesh.get(timestep, mesh.discretization_size() - 2),
                mesh.get(timestep - 1, mesh.discretization_size() - 1),
                k, flux));

            // After each run through, check that CFL satisfied.
            // TODO: remove redundant CFL checks.
            this->cfl_check_row(mesh, flux, delta_t, delta_x, timestep);
        }
    }

    /*
     * Stencils
     */
    static T leapfrog_stencil(T u_x_plus_1, T u_x_minus_1, T u_x_prev, double k, FluxFunction<T> *flux) {
        return u_x_prev - (flux->flux(u_x_plus_1) - flux->flux(u_x_minus_1)) * k;
    }
};

#endif //PDENCLOSE_LEAPFROGSOLVER_H