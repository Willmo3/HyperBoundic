//
// Created by will on 9/11/25.
//

#ifndef PDEAPPROX_LAXFRIEDRICHSSOLVER_H
#define PDEAPPROX_LAXFRIEDRICHSSOLVER_H
#include "PdeDiscretization.h"

/**
 * @brief Given a set of initial conditions over some discretization of a 1d space, a time discretization, and a number of timesteps,
 * Approximate the values of the system at different points in time.
 *
 * @param initial_state Array of size discretization_size representing the initial state of the system
 * @param num_timesteps Number of timesteps for the approximation.
 * @param delta_t Time discretization... i.e. how much time is passing logically for each step.
 * @param delta_x Space discretization... i.e. how much space is passing logically for each step.
 * @return a discretization of the partial differential equation system.
 */
PdeDiscretization lax_friedrichs_solver(const double* initial_state, uint32_t discretization_size, uint32_t num_timesteps, double delta_t, double delta_x);

#endif //PDEAPPROX_LAXFRIEDRICHSSOLVER_H