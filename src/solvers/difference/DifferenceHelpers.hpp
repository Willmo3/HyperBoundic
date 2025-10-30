//
// Created by will on 10/27/25.
//

#ifndef PDENCLOSE_DIFFERENCEHELPERS_H
#define PDENCLOSE_DIFFERENCEHELPERS_H
#include "domains/Numeric.hpp"
#include "meshes/RectangularMesh.hpp"
#include "meshes/CflCheck.hpp"

/**
 * @brief Perform the CFL check on a row of data.
 *
 * @param mesh Mesh to perform check on.
 * @param flux Flux function to use for CFL check.
 * @param delta_t Logical time discretization of the solver.
 * @param delta_x Logical space discretization of the solver.
 * @param timestep Timestep to evaluate.
 */
template<typename T>
requires Numeric<T>
static void cfl_check_row(const RectangularMesh<T> &mesh, FluxFunction<T> *flux, double delta_t, double delta_x, int timestep) {
    assert(timestep >= 0 && timestep < mesh.num_timesteps());
    for (auto point = 0; point < mesh.discretization_size(); point++) {
        if (!cfl_check(flux, mesh.get(timestep, point), delta_t, delta_x)) {
            std::cerr << "System blowup detected, exiting." << std::endl;
            exit(1);
        }
    }
}


#endif //PDENCLOSE_DIFFERENCEHELPERS_H