
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>

#include "solvers/LaxFriedrichsSolver.hpp"
#include "../Winterval/src/Winterval.h"
#include "domains/Real.hpp"

/**
 * Test a Lax Friedrichs approximation over the scalar and interval domains.
 */
void test_lf_unit_flux() {
    /*
     * Shared constants
     */
}

void test_lf_scalar() {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_x = 1; // Spatial discretization. Assuming total space of four split into 4 parts = 1.
    double delta_t = 1; // Spacing of time. Assuming operating over 4 logical time split into 4 timesteps = 1.

    auto initial_conditions = static_cast<Real *>(calloc(discretization_size, sizeof(double)));
    assert(initial_conditions);

    initial_conditions[0] = 1.0;
    initial_conditions[1] = 2.0;
    initial_conditions[2] = 3.0;
    initial_conditions[3] = 4.0;

    auto solution_matrix = LaxFriedrichsSolver<Real>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, LaxFriedrichsSolver<Real>::cubic_flux);
    solution_matrix.print_system();

    free(initial_conditions);
    initial_conditions = nullptr;
}

void test_lf_interval() {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_x = 1;
    double delta_t = 1;

    auto initial_conditions = static_cast<Winterval *>(calloc(discretization_size, sizeof(Winterval)));
    assert(initial_conditions);

    initial_conditions[0] = Winterval(0, 1);
    initial_conditions[1] = Winterval(1, 2);
    initial_conditions[2] = Winterval(2, 3);
    initial_conditions[3] = Winterval(3, 4);

    auto solution_matrix = LaxFriedrichsSolver<Winterval>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, LaxFriedrichsSolver<Winterval>::cubic_flux);
    solution_matrix.print_system();

    free(initial_conditions);
    initial_conditions = nullptr;
}

int main() {
    std::cout << "Scalar:" << std::endl;
    test_lf_scalar();
    std::cout << std::endl << "Interval:" << std::endl;
    test_lf_interval();
    return 0;
}