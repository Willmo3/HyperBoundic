//
// Created by will on 9/11/25.
//

#ifndef PDEAPPROX_SYSTEMAPPROXIMATION_H
#define PDEAPPROX_SYSTEMAPPROXIMATION_H
#include <array>
#include <cstdint>

/**
 * Discretization of a physical system represented by a partial differential equation.
 */
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
     */
    PdeDiscretization(uint32_t discretization_size, uint32_t num_timesteps, const double *initial_conditions);
    ~PdeDiscretization();

    /*
     * Accessors
     */
    uint32_t discretization_size() const;
    uint32_t num_timesteps() const;
    double get(uint32_t timestep, uint32_t index) const;
    void set(uint32_t timestep, uint32_t index, double value);

    /*
     * Helpers
     */
    void print_system() const;
private:
    double *system;
    const uint32_t _discretization_size;
    const uint32_t _num_timesteps;
};


#endif //PDEAPPROX_SYSTEMAPPROXIMATION_H