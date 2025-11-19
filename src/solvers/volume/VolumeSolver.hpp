//
// Created by will on 11/19/25.
//

#ifndef PDENCLOSE_VOLUMESOLVER_H
#define PDENCLOSE_VOLUMESOLVER_H
#include "domains/Numeric.hpp"
#include "flux/FluxFunction.hpp"
#include "meshes/RectangularMesh.hpp"

/**
 * Interface for finite volume method solver.
 * @tparam T Numeric type to operate over.
 */
template<typename T>
requires Numeric<T>
class VolumeSolver {
public:
    virtual ~VolumeSolver() = default;

    /**
     * @brief Approximate a finite volume mesh of a discretized system with a finite volume solver.
     * The mesh may be irregular, so the discretization constant is calculated for each control volume cell.
     *
     * @param initial_state Initial values of system.
     * @param width_values Width for each point in the initial state of the system -- irregular mesh.
     * This subsumes delta_x in the standard finite difference solver.
     * @param discretization_size Number of control volume cells.
     * @param num_timesteps Number of timesteps for simulation.
     * @param delta_t Change in time at each point.
     * @param flux Flux function to approximate system with.
     * @return The approximation of the system.
     */
    virtual RectangularMesh<T> solve(
        const std::vector<T> &initial_state,
        const std::vector<double> &width_values,
        uint32_t discretization_size,
        uint32_t num_timesteps,
        double delta_t,
        FluxFunction<T>* flux) = 0;

    /**
     * @brief Perform a CFL check over an entire solution mesh.
     * If fails, prints out the timestep and point of failure.
     *
     * @param solution Solution to check over.
     * @param flux Flux function to check CFL satisfiability.
     * @param delta_t Temporal discretization constant.
     * @param width_values Spatial discretization values -- different for each mesh point.
     *
     * @return Whether the CFL check passed for the entire mesh.
     */
    bool cfl_check_mesh(const RectangularMesh<T> &solution, FluxFunction<T> *flux, double delta_t, const std::vector<double> &width_values) {
        assert(width_values.size() == solution.discretization_size());

        for (auto timestep = 0; timestep < solution.num_timesteps(); timestep++) {
            for (auto point = 0; point < solution.discretization_size(); point++) {
                if (!cfl_check(flux, solution.get(timestep, point), delta_t, width_values[point])) {
                    std::cout << "First CFL violation at timestep " << timestep << ", point " << point << std::endl;
                    return false;
                }
            }
        }
        std::cout << "No CFL violations found." << std::endl;
        return true;
    }

};

#endif //PDENCLOSE_VOLUMESOLVER_H