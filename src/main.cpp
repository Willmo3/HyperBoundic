
#include <fstream>

#include "meshes/RectangularMesh.hpp"
#include "domains/Real.hpp"
#include "solvers/flux/CubicFlux.hpp"
#include "solvers/volume/LocalLaxFriedrichsSolver.hpp"
#include "Wixed/WixedForm.hpp"

void test_llf_real() {
    auto discretization_size = 4;
    auto num_timesteps = 4;

    auto initial_conditions = std::unique_ptr<Real>(static_cast<Real *>(calloc(discretization_size, sizeof(double))));
    initial_conditions.get()[0] = 1.0;
    initial_conditions.get()[1] = 2.0;
    initial_conditions.get()[2] = 3.0;
    initial_conditions.get()[3] = 4.0;

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

int main() {
    test_llf_real();

    return 0;
}
