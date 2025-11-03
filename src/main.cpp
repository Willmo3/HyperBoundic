
#include "meshes/RectangularMesh.hpp"
#include "domains/Real.hpp"
#include "solvers/flux/CubicFlux.hpp"
#include "solvers/volume/LocalLaxFriedrichsSolver.hpp"
#include "DualDomain/MixedForm.hpp"
#include "io/save_initial_conditions.hpp"
#include "solvers/difference/LaxFriedrichsSolver.hpp"
#include "solvers/flux/BuckleyLeverettFlux.hpp"
#include "solvers/flux/BurgersFlux.hpp"
#include "solvers/flux/LwrFlux.hpp"

void test_llf_real() {
    auto discretization_size = 4;
    auto num_timesteps = 4;

    auto initial_conditions = std::vector<Real>(discretization_size);
    initial_conditions[0] = 1.0;
    initial_conditions[1] = 2.0;
    initial_conditions[2] = 3.0;
    initial_conditions[3] = 4.0;

    // Mesh more fine grained in the center.
    auto delta_t = 0.02;
    auto width_values = std::vector<double>(discretization_size);
    width_values[0] = 1;
    width_values[1] = 1;
    width_values[2] = 1;
    width_values[3] = 1;

    // TODO: free flux fns when finished.
    auto solution_matrix = LocalLaxFriedrichsSolver<Real>::solve(initial_conditions, width_values, discretization_size, num_timesteps, delta_t, new CubicFlux<Real>);
    solution_matrix.print_system();
}

void test_flux() {
    uint32_t discretization_size = 4;
    uint32_t num_timesteps = 30;
    double delta_t = 0.02; // Spacing of time. Assuming operating over 4 logical time split into 4 timesteps = 1.
    double delta_x = 1; // Spatial discretization. Assuming total space of four split into 4 parts = 1.

    auto initial_conditions = std::vector<Real>(discretization_size);

    initial_conditions[0] = 1;
    initial_conditions[1] = 2;
    initial_conditions[2] = 3;
    initial_conditions[3] = 4;

    auto solution_matrix = LaxFriedrichsSolver<Real>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new BuckleyLeverett<Real>());
    solution_matrix.print_system();
}

int main() {
    std::vector<Real> initial_conditions = std::vector<Real>(4);
    initial_conditions[0] = 1;
    initial_conditions[1] = 2;
    initial_conditions[2] = 3;
    initial_conditions[3] = 4;

    write_initial_conditions("initial_conditions.json", initial_conditions);
    return 0;
}
