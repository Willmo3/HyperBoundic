
#include "meshes/RectangularMesh.hpp"
#include "domains/Real.hpp"
#include "flux/CubicFlux.hpp"
#include "solvers/volume/LocalLaxFriedrichsSolver.hpp"
#include "DualDomain/MixedForm.hpp"
#include "io/initial_conditions.hpp"
#include "io/SimulationConfig.hpp"
#include "solvers/difference/LaxFriedrichsSolver.hpp"
#include "flux/BuckleyLeverettFlux.hpp"
#include "flux/BurgersFlux.hpp"
#include "flux/LwrFlux.hpp"
#include "args/match_names.hpp"
#include "io/generate_source_files.h"

/**
 * Run a small simulation with 4 timesteps.
 * @param cfg_path Path to configuration file
 * @param initial_conds_path Path to string with initial conditions.
 */
void run_simulation(const std::string &cfg_path, const std::string &initial_conds_path);
/**
 *
 * @param argc Number of arguments
 * @param argv Argument vector
 * @param write_test Pointer to option about whether to write out a sanity test file.
 * @param cfg_path Pointer to string where path of discretization config will be placed.
 * @param initial_conds_path Pointer to string where path of initial conditions will be placed.
 * @return whether no invalid arguments were provided
 */
static bool get_args(int argc, char *argv[], bool *write_test, std::string *cfg_path, std::string *initial_conds_path);

/**
 * Print usage information to stdout.
 */
static void usage();

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

// Note: for now, assume only real-valued.
int main(int argc, char *argv[]) {
    std::string cfg_path = "";
    std::string initial_conds_path = "";
    bool gen_sources = false;

    if (argc == 1) {
        std::cout << "No arguments provided, running sanity test." << std::endl;
        // replace with sanity test
        test_llf_real();
        exit(EXIT_SUCCESS);
    }

    // Read command line args.
    if (!get_args(argc, argv, &gen_sources, &cfg_path, &initial_conds_path)) {
        std::cerr << "Invalid arguments." << std::endl;
        usage();
        exit(EXIT_FAILURE);
    }

    // Validate command line args.
    if (gen_sources != cfg_path.empty()) {
        std::cerr << "Either write a sanity test or run a simulation." << std::endl;
        usage();
        exit(EXIT_FAILURE);
    }
    if (cfg_path.empty() != initial_conds_path.empty()) {
        std::cerr << "Specify both an initial conditions file and a configuration file." << std::endl;
        usage();
        exit(EXIT_FAILURE);
    }

    if (gen_sources) {
        generate_source_files();
    } else {
        run_simulation(cfg_path, initial_conds_path);
    }

    return 0;
}

static void usage() {
    std::cout << R"(Usage: "PDEnclose -w" OR "PDEnclose -c <config_path> -s <initial_conditions_path>")" << std::endl;
    std::cout << "\t-w: Write out source files for testing." << std::endl;
    std::cout << "\t-c: Path to configuration file." << std::endl;
    std::cout << "\t-s: Path to initial conditions file." << std::endl;
}

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

void run_simulation(const std::string &cfg_path, const std::string &initial_conds_path) {
    // Read config
    auto config = read_config(cfg_path);

    // Note: duplicated logic between these branches because we can only initialize once we know the templated type.
    if (config.domain == "real") {
        auto initial_conditions = read_initial_conditions<Real>(initial_conds_path);
        auto flux = match_flux<Real>(config.flux);
        // For now, only difference solvers.
        auto solver = match_difference<Real>(config.solver);
        auto solution = solver->solve(initial_conditions, config.discretization_size, config.num_timesteps, config.delta_t, config.delta_x, flux);
        solution.print_system();

        delete solver;
        delete flux;

    } else if (config.domain == "interval") {
        auto initial_conditions = read_initial_conditions<Winterval>(initial_conds_path);
        auto flux = match_flux<Winterval>(config.flux);
        // For now, only difference solvers.
        auto solver = match_difference<Winterval>(config.solver);
        auto solution = solver->solve(initial_conditions, config.discretization_size, config.num_timesteps, config.delta_t, config.delta_x, flux);
        solution.print_system();

        delete solver;
        delete flux;

    } else if (config.domain == "affine") {
        auto initial_conditions = read_initial_conditions<AffineForm>(initial_conds_path);
        auto flux = match_flux<AffineForm>(config.flux);
        // For now, only difference solvers.
        auto solver = match_difference<AffineForm>(config.solver);
        auto solution = solver->solve(initial_conditions, config.discretization_size, config.num_timesteps, config.delta_t, config.delta_x, flux);
        solution.print_system();

        delete solver;
        delete flux;

    } else if (config.domain == "mixed") {
        auto initial_conditions = read_initial_conditions<MixedForm>(initial_conds_path);
        auto flux = match_flux<MixedForm>(config.flux);
        // For now, only difference solvers.
        auto solver = match_difference<MixedForm>(config.solver);
        auto solution = solver->solve(initial_conditions, config.discretization_size, config.num_timesteps, config.delta_t, config.delta_x, flux);
        solution.print_system();

        delete solver;
        delete flux;

    } else {
        std::cerr << "Invalid domain!" << std::endl;
        exit(EXIT_FAILURE);
    }
}
