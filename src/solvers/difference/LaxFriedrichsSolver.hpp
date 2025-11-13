//
// Created by will on 9/11/25.
//

#ifndef PDENCLOSE_LAXFRIEDRICHSSOLVER_H
#define PDENCLOSE_LAXFRIEDRICHSSOLVER_H

#include <cmath>

#include "domains/Numeric.hpp"
#include "meshes/Rectangularmesh.hpp"
#include "flux/FluxFunction.hpp"
#include "DifferenceSolver.hpp"

template<typename T>
requires Numeric<T>
class LaxFriedrichsSolver final : public DifferenceSolver<T> {

public:
    /*
     * Constructor
     */
    LaxFriedrichsSolver(const std::vector<T> &initial_state,
                        uint32_t discretization_size,
                        uint32_t num_timesteps,
                        double delta_t,
                        double delta_x,
                        FluxFunction<T> *flux):
        DifferenceSolver<T>(initial_state, discretization_size, num_timesteps, delta_t, delta_x, flux) {}

    /*
     * Super fields
     */
    using DifferenceSolver<T>::delta_t;
    using DifferenceSolver<T>::delta_x;
    using DifferenceSolver<T>::mesh;
    using DifferenceSolver<T>::flux;

    /*
     * Solver
     */
    void solve() override {
        auto k = delta_t / delta_x * 1/2;

        for (auto timestep = 0; timestep < mesh.num_timesteps() - 1; timestep++) {
            for (auto x = 1; x < mesh.discretization_size() - 1; x++) {
                auto u_x_plus_1 = mesh.get(timestep, x + 1);
                auto u_x_minus_1 = mesh.get(timestep, x - 1);
                mesh.set(timestep + 1, x, lax_friedrichs_stencil(u_x_plus_1, u_x_minus_1, k, flux));
            }

            // Currently, only support periodic boundary conditions.
            mesh.set(timestep + 1, 0, lax_friedrichs_stencil(
                mesh.get(timestep, 1),
                mesh.get(timestep, mesh.discretization_size() - 1),
                k, flux));

            mesh.set(timestep + 1, mesh.discretization_size() - 1, lax_friedrichs_stencil(
                mesh.get(timestep, 0),
                mesh.get(timestep, mesh.discretization_size() - 2),
                k, flux));

            // After each run through, check that CFL satisfied.
            // TODO: remove redundant CFL checks.
            this->cfl_check_row(mesh, flux, delta_t, delta_x, timestep);
        }
    }

    /*
     * Stencils
     */
    static T lax_friedrichs_stencil(T u_i_plus_1, T u_i_minus_1, double k, FluxFunction<T> *flux) {
        return (u_i_plus_1 + u_i_minus_1) * 0.5 - (flux->flux(u_i_plus_1) - flux->flux(u_i_minus_1)) * k;
    }
};

#endif //PDENCLOSE_LAXFRIEDRICHSSOLVER_H