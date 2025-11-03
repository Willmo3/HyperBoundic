
#include "../meshes/RectangularMesh.hpp"
#include "../domains/Real.hpp"
#include "../solvers/flux/CubicFlux.hpp"
#include "../solvers/volume/LocalLaxFriedrichsSolver.hpp"
#include "DualDomain/MixedForm.hpp"
#include "initial_conditions.hpp"
#include "SimulationConfig.hpp"
#include "../solvers/difference/LaxFriedrichsSolver.hpp"
#include "../solvers/flux/BuckleyLeverettFlux.hpp"
#include "../solvers/flux/BurgersFlux.hpp"
#include "../solvers/flux/LwrFlux.hpp"

// void test_llf_real() {
//     auto discretization_size = 4;
//     auto num_timesteps = 4;
//
//     auto initial_conditions = std::vector<Real>(discretization_size);
//     initial_conditions[0] = 1.0;
//     initial_conditions[1] = 2.0;
//     initial_conditions[2] = 3.0;
//     initial_conditions[3] = 4.0;
//
//     // Mesh more fine grained in the center.
//     auto delta_t = 0.02;
//     auto width_values = std::vector<double>(discretization_size);
//     width_values[0] = 1;
//     width_values[1] = 1;
//     width_values[2] = 1;
//     width_values[3] = 1;
//
//     // TODO: free flux fns when finished.
//     auto solution_matrix = LocalLaxFriedrichsSolver<Real>::solve(initial_conditions, width_values, discretization_size, num_timesteps, delta_t, new CubicFlux<Real>);
//     solution_matrix.print_system();
// }
//
// void test_flux() {
//     uint32_t discretization_size = 4;
//     uint32_t num_timesteps = 30;
//     double delta_t = 0.02; // Spacing of time. Assuming operating over 4 logical time split into 4 timesteps = 1.
//     double delta_x = 1; // Spatial discretization. Assuming total space of four split into 4 parts = 1.
//
//     auto initial_conditions = std::vector<Real>(discretization_size);
//
//     initial_conditions[0] = 1;
//     initial_conditions[1] = 2;
//     initial_conditions[2] = 3;
//     initial_conditions[3] = 4;
//
//     auto solution_matrix = LaxFriedrichsSolver<Real>::solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new BuckleyLeverett<Real>());
//     solution_matrix.print_system();
// }


/**
 * Write a small sanity test of initial conditions.
 */
void write_sanity_conditions() {
    // Write input config.
    auto cfg = SimulationConfig("burgers", "real", 4, 4, 1, 0.01);

    // Write actual initial conditions.
    std::vector<Real> initial_conditions = std::vector<Real>(4);
    initial_conditions[0] = 1;
    initial_conditions[1] = 2;
    initial_conditions[2] = 3;
    initial_conditions[3] = 4;

    write_initial_conditions("sanity_conditions.json", initial_conditions);
    write_config("sanity_config.json", cfg);
}

/**
 * Run a small simulation with 4 timesteps.
 * @param cfg_path Path to configuration file.
 * @param initial_conds_path Path to string with initial conditions.
 */
void run_simulation(const std::string &cfg_path, const std::string &initial_conds_path) {
    // Read config
    auto config = read_config(cfg_path);

    // TODO: dictionary mapping strings to fluxes.

    // TODO: generate flux function
    if (config.domain == "real") {
        auto initial_conditions = read_initial_conditions<Real>(initial_conds_path);
        // TODO: require delta_t, delta_x
        auto solution = LaxFriedrichsSolver<Real>::solve(initial_conditions, config.discretization_size, config.num_timesteps, config.delta_t, config.delta_x, new CubicFlux<Real>());
        solution.print_system();
    } else {
        std::cerr << "Invalid domain!" << std::endl;
        exit(EXIT_FAILURE);
    }
}

/**
 *
 * @param argc Number of arguments
 * @param argv Argument vector
 * @param write_test Pointer to option about whether to write out a sanity test file.
 * @param cfg_path Pointer to string where path of discretization config will be placed.
 * @param initial_conds_path Pointer to string where path of initial conditions will be placed.
 * @return whether no invalid arguments were provided
 */
static bool get_args(int argc, char *argv[], bool *write_test, std::string *cfg_path, std::string *initial_conds_path) {
    int ch = 0;
    while ((ch = getopt(argc, argv, "wc:s:")) != -1) {
        switch (ch) {
            case 'w':
                *write_test = true;
                break;
            case 's':
                *initial_conds_path = optarg;
                break;
            case 'c':
                *cfg_path = optarg;
                break;
            default:
                return false;
        }
    }
    return true;
}

// Note: for now, assume only real-valued.
int main(int argc, char *argv[]) {
    std::string cfg_path = "";
    std::string initial_conds_path = "";
    bool write_test = false;

    // Read command line args.
    if (!get_args(argc, argv, &write_test, &cfg_path, &initial_conds_path)) {
        std::cerr << "Invalid arguments." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Validate command line args.
    if (write_test != cfg_path.empty()) {
        std::cerr << "Either write a sanity test or run a simulation." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (cfg_path.empty() != initial_conds_path.empty()) {
        std::cerr << "Specify both an initial conditions file and a configuration file." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (write_test) {
        // Put out a small sanity test.
        write_sanity_conditions();
    } else {
        run_simulation(cfg_path, initial_conds_path);
    }

    return 0;
}
