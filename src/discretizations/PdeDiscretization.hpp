//
// Created by will on 9/11/25.
//

#ifndef PDEAPPROX_SYSTEMAPPROXIMATION_H
#define PDEAPPROX_SYSTEMAPPROXIMATION_H
#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>

#include "../solvers/flux/FluxFunction.hpp"
#include "../domains/Numeric.hpp"

/**
 * Discretization of a physical system represented by a partial differential equation.
 * @param T numeric type to approximate system.
 */
template<typename T>
requires Numeric<T>
class PdeDiscretization {
public:
    /*
     * Constructors
     */

    /**
     *
     * @param discretization_size Size of discretization, > 0.
     * @param num_timesteps Number of timesteps for this discretization.
     * @param initial_conditions Array of starting conditions for the system, of len discretization_size.
     * This must be a unique pointer, since values will be consumed by this function.
     */
    PdeDiscretization(uint32_t discretization_size, uint32_t num_timesteps, std::unique_ptr<T> &initial_conditions)
        :_discretization_size(discretization_size),  _num_timesteps(num_timesteps){
        assert(discretization_size > 0);
        assert(num_timesteps > 0);

        system = static_cast<T *>(calloc(sizeof(T), num_timesteps * discretization_size));
        assert(system);

        // Copy initial values into solution matrix
        // TODO: should be doable with memcpy
        for (auto i = 0; i < discretization_size; i++) {
            system[i] = initial_conditions.get()[i];
        }
    }
    ~PdeDiscretization() {
        free(system);
        system = nullptr;
    }

    /*
     * Accessors
     */
    uint32_t discretization_size() const {
        return _discretization_size;
    }
    uint32_t num_timesteps() const {
        return _num_timesteps;
    }
    T get(uint32_t timestep, uint32_t index) const {
        assert(timestep < _num_timesteps);
        assert(index < _discretization_size);
        return system[timestep * _discretization_size + index];
    }
    void set(uint32_t timestep, uint32_t index, T value) {
        assert(timestep < _num_timesteps);
        assert(index < _discretization_size);
        system[timestep * _discretization_size + index] = value;
    }

    /*
     * Validators
     */

    /**
     *
     * @param flux Flux function to use for CFL check.
     * @param delta_t Logical time discretization of the solver.
     * @param delta_x Logical space discretization of the solver.
     * @param timestep Timestep to evaluate.
     */
    void cfl_check_row(FluxFunction<T> *flux, double delta_t, double delta_x, int timestep) {
        assert(timestep >= 0 && timestep < _num_timesteps);
        for (auto point = 0; point < _discretization_size; point++) {
            if (!cfl_check(flux, get(timestep, point), delta_t, delta_x)) {
                std::cerr << "System blowup detected, exiting." << std::endl;
                exit(1);
            }
        }
    }

    /*
     * Helpers
     */
    void print_system() const {
        for (auto t = 0; t < _num_timesteps; t++) {
            std::cout << "T" << t << ": ";
            for (auto i = 0; i < _discretization_size; i++) {
                std::cout << system[_discretization_size * t + i] << " ";
            }
            std::cout << std::endl;
        }
    }
private:
    T *system;
    const uint32_t _discretization_size;
    const uint32_t _num_timesteps;
};

#endif //PDEAPPROX_SYSTEMAPPROXIMATION_H