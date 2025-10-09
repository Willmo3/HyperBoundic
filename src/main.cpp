
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>

#include "solvers/LaxFriedrichsSolver.hpp"
#include "../lib/Winterval/Winterval.hpp"
#include "../lib/Waffine/WaffineForm.hpp"
#include "../lib/Wixed/WixedForm.hpp"
#include "domains/Real.hpp"
#include "solvers/flux/CubicFlux.hpp"
#include "visualization/DiscretizationVisualizers.hpp"

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
    double delta_t = 0.02; // Spacing of time. Assuming operating over 4 logical time split into 4 timesteps = 1.

    auto initial_conditions = std::unique_ptr<Real>(static_cast<Real *>(calloc(discretization_size, sizeof(double))));
    assert(initial_conditions);

    initial_conditions.get()[0] = 1.0;
    initial_conditions.get()[1] = 2.0;
    initial_conditions.get()[2] = 3.0;
    initial_conditions.get()[3] = 4.0;

    auto solution_matrix = LaxFriedrichsSolver<Real>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<Real>());
    solution_matrix.print_system();
    // show_real_surface(&solution_matrix);
}

void test_lf_interval() {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_x = 1;
    double delta_t = 0.02;

    auto initial_conditions = std::unique_ptr<Winterval>(static_cast<Winterval *>(calloc(discretization_size, sizeof(Winterval))));

    initial_conditions.get()[0] = Winterval(0, 1);
    initial_conditions.get()[1] = Winterval(1, 2);
    initial_conditions.get()[2] = Winterval(2, 3);
    initial_conditions.get()[3] = Winterval(3, 4);

    auto solution_matrix = LaxFriedrichsSolver<Winterval>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<Winterval>());
    solution_matrix.print_system();
}
// TODO: leapfrog, reduced product

void test_lf_affine() {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_x = 1;
    double delta_t = 0.02;

    auto initial_conditions = std::unique_ptr<WaffineForm>(static_cast<WaffineForm *>(calloc(discretization_size, sizeof(WaffineForm))));

    initial_conditions.get()[0] = WaffineForm(Winterval(0, 1));
    initial_conditions.get()[1] = WaffineForm(Winterval(1, 2));
    initial_conditions.get()[2] = WaffineForm(Winterval(2, 3));
    initial_conditions.get()[3] = WaffineForm(Winterval(3, 4));

    auto solution_matrix = LaxFriedrichsSolver<WaffineForm>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<WaffineForm>());
    solution_matrix.print_system();
}

void test_lf_mixed() {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 4;
    double delta_x = 1;
    double delta_t = 0.02;

    auto initial_conditions = std::unique_ptr<WixedForm>(static_cast<WixedForm *>(calloc(discretization_size, sizeof(WixedForm))));

    initial_conditions.get()[0] = WixedForm(Winterval(0, 1));
    initial_conditions.get()[1] = WixedForm(Winterval(1, 2));
    initial_conditions.get()[2] = WixedForm(Winterval(2, 3));
    initial_conditions.get()[3] = WixedForm(Winterval(3, 4));

    auto solution_matrix = LaxFriedrichsSolver<WixedForm>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new CubicFlux<WixedForm>());
    solution_matrix.print_system();
}

int main() {
    std::cout << "Scalar:" << std::endl;
    test_lf_scalar();
    std::cout << std::endl << "Interval:" << std::endl;
    test_lf_interval();
    std::cout << std::endl << "Affine:" << std::endl;
    test_lf_affine();
    std::cout << std::endl << "Mixed:" << std::endl;
    test_lf_mixed();
    return 0;
}