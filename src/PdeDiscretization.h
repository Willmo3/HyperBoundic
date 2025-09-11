//
// Created by will on 9/11/25.
//

#ifndef PDEAPPROX_SYSTEMAPPROXIMATION_H
#define PDEAPPROX_SYSTEMAPPROXIMATION_H
#include <array>
#include <cstdint>


class PdeDiscretization {
public:
    /*
     * Constructors
     */
    PdeDiscretization(uint32_t discretization_size, uint32_t num_timesteps, const double *constinitial_conditions);
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