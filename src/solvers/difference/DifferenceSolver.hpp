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
    /**
     * @param initial_state vector of size discretization_size representing the initial state of the system
     * @param num_timesteps Number of timesteps for the approximation. (num rows)
     * @param discretization_size Size of space being approximated. (num cols)
     * @param delta_t Time discretization... i.e. how much time is passing logically for each step. Must be > 0, < INFINITY
     * @param delta_x Space discretization... i.e. how much space is passing logically for each step. Must be > 0, < INFINITY
     * @param flux Flux function to use for this approximation.
     */
    DifferenceSolver(const std::vector<T> &initial_state, uint32_t discretization_size, uint32_t num_timesteps, double delta_t, double delta_x, FluxFunction<T> *flux): // Copy initial state
        flux(flux), delta_t(delta_t), delta_x(delta_x), mesh(RectangularMesh<T>(discretization_size, num_timesteps)) {

        assert(delta_t > 0 && delta_t < INFINITY);
        assert(delta_x > 0 && delta_x < INFINITY);
        mesh.copy_initial_conditions(initial_state);
    }

    DifferenceSolver() = default;
    virtual ~DifferenceSolver() {
        // delete flux;
    }

    /**
     * @brief Given a set of initial conditions over some discretization of a 1d space, a time discretization, and a number of timesteps,
     * Approximate the values of the system at different points in time.
     */
    virtual void solve() = 0;

    const RectangularMesh<T> &get_mesh() {
        return mesh;
    }

protected:
    // TODO: Remove unneeded params
    void cfl_check_row(const RectangularMesh<T> &mesh, FluxFunction<T> *flux, double delta_t, double delta_x, int timestep) {
        assert(timestep >= 0 && timestep < mesh.num_timesteps());
        for (auto point = 0; point < mesh.discretization_size(); point++) {
            if (!cfl_check(flux, mesh.get(timestep, point), delta_t, delta_x)) {
                std::cerr << "System blowup detected at timestep " << timestep << ", point " << point << std::endl;
            }
        }
    }

    RectangularMesh<T> mesh;
    FluxFunction<T> *flux;
    double delta_t;
    double delta_x;
};

#endif //PDENCLOSE_PDESOLVER_H
