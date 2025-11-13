//
// Created by will on 11/4/25.
//

#ifndef PDENCLOSE_PDESOLVER_H
#define PDENCLOSE_PDESOLVER_H
#include "domains/Numeric.hpp"
#include "flux/FluxFunction.hpp"
#include "meshes/RectangularMesh.hpp"
#include "meshes/CflCheck.hpp"

/**
 * Interface for finite difference method solver.
 * @tparam T Numeric type being solved over.
 */
template<typename T>
requires Numeric<T>
class DifferenceSolver {
public:
    virtual ~DifferenceSolver() = default;

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
    virtual RectangularMesh<T> solve(
        const std::vector<T> &initial_state,
        uint32_t discretization_size,
        uint32_t num_timesteps,
        double delta_t,
        double delta_x,
        FluxFunction<T>* flux) = 0;

    /**
     * @brief Perform a CFL check over an entire solution mesh.
     *
     * @param solution Solution to check over.
     * @param flux Flux function to check CFL satisfiability.
     * @param delta_t Temporal discretization constant.
     * @param delta_x Spatial discretization constant.
     */
    void    cfl_check_mesh(const RectangularMesh<T> &solution, FluxFunction<T> *flux, double delta_t, double delta_x) {
        for (auto timestep = 0; timestep < solution.num_timesteps(); timestep++) {
            for (auto point = 0; point < solution.discretization_size(); point++) {
                if (!cfl_check(flux, solution.get(timestep, point), delta_t, delta_x)) {
                    std::cerr << "System blowup detected at timestep " << timestep << ", point " << point << std::endl;
                }
            }
        }
    }
};

#endif //PDENCLOSE_PDESOLVER_H
