//
// Created by will on 9/11/25.
//

#include "PdeDiscretization.h"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <ostream>

/*
 * Constructors
 */
PdeDiscretization::PdeDiscretization(uint32_t discretization_size, uint32_t num_timesteps, const double* const initial_conditions)
    :_discretization_size(discretization_size),  _num_timesteps(num_timesteps) {

    system = static_cast<double *>(calloc(sizeof(double), num_timesteps * discretization_size));
    assert(system);

    // Copy initial values into solution matrix
    assert(memcpy(system, initial_conditions, discretization_size * sizeof(double)));
}
PdeDiscretization::~PdeDiscretization() {
    free(system);
    system = nullptr;
}

/*
 * Accessors
 */
uint32_t PdeDiscretization::discretization_size() const {
    return _discretization_size;
}
uint32_t PdeDiscretization::num_timesteps() const {
    return _num_timesteps;
}

double PdeDiscretization::get(uint32_t timestep, uint32_t index) const {
    assert(timestep < _num_timesteps);
    assert(index < _discretization_size);
    return system[timestep * _discretization_size + index];
}
void PdeDiscretization::set(uint32_t timestep, uint32_t index, double value) {
    assert(timestep < _num_timesteps);
    assert(index < _discretization_size);
    system[timestep * _discretization_size + index] = value;
}

/*
 * Helpers
 */
void PdeDiscretization::print_system() const {
    for (auto t = 0; t < _num_timesteps; t++) {
        std::cout << "T" << t << ": ";
        for (auto i = 0; i < _discretization_size; i++) {
            std::cout << system[_discretization_size * t + i] << " ";
        }
        std::cout << std::endl;
    }
}



