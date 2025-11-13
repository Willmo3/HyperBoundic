//
// Created by will on 10/27/25.
//

#ifndef PDENCLOSE_LOCALLAXFRIEDRICHSSOLVER_H
#define PDENCLOSE_LOCALLAXFRIEDRICHSSOLVER_H
#include <cmath>

#include "meshes/RectangularMesh.hpp"
#include "meshes/CflCheck.hpp"

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
                auto u_x_plus_1 = solution.get(timestep, x + 1);
                auto u_x_minus_1 = solution.get(timestep, x - 1);

                auto k = viscosity_coefficient(u_x_plus_1, u_x_minus_1, flux) * 1/2;
                solution.set(timestep + 1, x, local_lax_friedrichs_stencil(u_x_plus_1, u_x_minus_1, k, flux));
            }

            // Currently, only support periodic boundary conditions.
            auto k_0 = viscosity_coefficient(solution.get(timestep, 1), solution.get(timestep, solution.discretization_size() - 1), flux) * 1/2;
            solution.set(timestep + 1, 0, local_lax_friedrichs_stencil(
                solution.get(timestep, 1),
                solution.get(timestep, solution.discretization_size() - 1),
                k_0, flux));

            auto k_n = viscosity_coefficient(solution.get(timestep, 0), solution.get(timestep, discretization_size - 2), flux) * 1/2;
            solution.set(timestep + 1, solution.discretization_size() - 1, local_lax_friedrichs_stencil(
                solution.get(timestep, 0),
                solution.get(timestep, discretization_size - 2),
                k_n, flux));

            // After each run through, check that CFL satisfied.
            // TODO: cfl_check_row_irregular?
            for (auto x = 0; x < solution.discretization_size(); x++) {
                if (!cfl_check(flux, solution.get(timestep, x), delta_t, width_values[x])) {;
                    std::cerr << "CFL check failed at point " << x << " at timestep " << timestep << "." << std::endl;
                }
            }
        }

        // CFL check last row
        for (auto x = 0; x < solution.discretization_size(); x++) {
            if (!cfl_check(flux, solution.get(num_timesteps - 1, x), delta_t, width_values[x])) {;
                std::cerr << "CFL check failed at point " << x << " at timestep " << (num_timesteps - 1) << "." << std::endl;
            }
        }

        return solution;
    }

    /*
     * In general, the viscosity of a cel is defined by the eigenvalues of the flux's Jacobian at the left and right states.
     * However, since we have a 1D system, this reduces to the absolute values of the derivatives at the left and right states.
     */
    static T viscosity_coefficient(T u_i_plus_1, T u_i_minus_1, FluxFunction<T> *flux) {
        auto right_propagation = flux->derivative_flux(u_i_plus_1).abs();
        auto left_propagation = flux->derivative_flux(u_i_minus_1).abs();
        return std::max(right_propagation, left_propagation);
    }

    // LLF stencil derived from: https://epubs.siam.org/doi/epdf/10.1137/0909030
    // (see eqn 4.11. Characteristic speeds in 1d care only about left and right)
    // The application in the 1d case is clearer in https://www.martin-schreiber.info/data/webdata/phd_thesis_html/schreiber14dissertationse12.html
    // See section 2.10.1 for example with Jacobians more clearly marked. Since they consider 2d, we can replace 1d case with scalar derivative.
    // Rusanov
    static T local_lax_friedrichs_stencil(T u_i_plus_1, T u_i_minus_1, T k, FluxFunction<T> *flux) {
        return (flux->flux(u_i_plus_1) + flux->flux(u_i_minus_1)) * 0.5 - (u_i_plus_1 - u_i_minus_1) * k * 0.5;
    }
};

#endif //PDENCLOSE_LOCALLAXFRIEDRICHSSOLVER_H