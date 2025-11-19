
#include "meshes/RectangularMesh.hpp"
#include "domains/Real.hpp"
#include "solvers/volume/LocalLaxFriedrichsSolver.hpp"
#include "DualDomain/MixedForm.hpp"
#include "experiment/SimulationConfig.hpp"
#include "solvers/difference/LaxFriedrichsSolver.hpp"
#include "flux/BurgersFlux.hpp"
#include "args/match_names.hpp"
#include "experiment/generators/generate_initial_conditions.hpp"
#include "experiment/generators/generate_source_files.h"
#include "visualization/MeshVisualizer.hpp"

/**
 * Run a user-configured simulation
 * @param cfg_path Path to configuration file
 * @param initial_conds_path Path to string with initial conditions.
 * @param run_cfl whether to run a CFL check pass on the solution mesh after execution.
 */
void run_simulation(const std::string &cfg_path, const std::string &initial_conds_path, bool run_cfl);
/**
 *
 * @param argc Number of arguments
 * @param argv Argument vector
 * @param write_test Pointer to option about whether to write out a sanity test file.
 * @param run_cfl Pointer to option about whether to run CFL check after executing simulation.
 * @param cfg_path Pointer to string where path of discretization config will be placed.
 * @param initial_conds_path Pointer to string where path of initial conditions will be placed.
 * @return whether no invalid arguments were provided
 */
static bool get_args(int argc, char *argv[], bool *write_test, bool *run_cfl, std::string *cfg_path, std::string *initial_conds_path);

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

    auto delta_t = 0.02;
    auto width_values = std::vector<double>(discretization_size);
    width_values[0] = 100;
    width_values[1] = 100;
    width_values[2] = 100;
    width_values[3] = 100;

    // TODO: free flux fns when finished.
    auto solution_matrix = LocalLaxFriedrichsSolver<Real>().solve(initial_conditions, width_values, discretization_size, num_timesteps, delta_t, new BurgersFlux<Real>);
    solution_matrix.print_system();
}

void show_pathological_burgers() {
    auto delta_x = 1.0;
    auto delta_t = 10.0;

    auto discretization_size = 5;
    auto num_timesteps = 5;

    auto initial_conditions = std::vector<Real>(discretization_size);
    initial_conditions[0] = 0.39;
    initial_conditions[1] = 0.66;
    initial_conditions[2] = 0.84;
    initial_conditions[3] = 0.75;
    initial_conditions[4] = 0.21;

    auto solution_matrix = LeapfrogSolver<Real>().solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new BurgersFlux<Real>());
    save_real_surface(&solution_matrix, "pathological_burgers_leapfrog.svg");
    // show_real_surface(&solution_matrix);
}

void show_acceptable_burgers() {
    auto delta_x = 1;
    auto delta_t = 1;

    auto discretization_size = 5;
    auto num_timesteps = 5;

    auto initial_conditions = std::vector<Real>(discretization_size);
    initial_conditions[0] = 0.39;
    initial_conditions[1] = 0.66;
    initial_conditions[2] = 0.84;
    initial_conditions[3] = 0.75;
    initial_conditions[4] = 0.21;

    auto solution_matrix = LeapfrogSolver<Real>().solve(initial_conditions, discretization_size, num_timesteps, delta_t, delta_x, new BurgersFlux<Real>());
    save_real_surface(&solution_matrix, "acceptable_burgers_leapfrog.svg");
    // show_real_surface(&solution_matrix);
}

// Note: for now, assume only real-valued.
int main(int argc, char *argv[]) {
    std::string cfg_path = "";
    std::string initial_conds_path = "";
    bool gen_sources = false;
    bool run_cfl = false;

    if (argc == 1) {
        std::cout << "No arguments provided, running sanity test." << std::endl;
        // replace with sanity test
        // test_llf_real();
        show_pathological_burgers();
        show_acceptable_burgers();
        exit(EXIT_SUCCESS);
    }

    // Read command line args.
    if (!get_args(argc, argv, &gen_sources, &run_cfl, &cfg_path, &initial_conds_path)) {
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
    if (gen_sources && run_cfl) {
        std::cerr << "Cannot run CFL check when generating source files." << std::endl;
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
        run_simulation(cfg_path, initial_conds_path, run_cfl);
    }

    return 0;
}

static void usage() {
    std::cout << R"(Usage: "PDEnclose -w" OR "PDEnclose -c <config_path> -s <initial_conditions_path> [-t]")" << std::endl;
    std::cout << "\t-w: Write out source files for testing." << std::endl;
    std::cout << "\t-c: Path to configuration file." << std::endl;
    std::cout << "\t-s: Path to initial conditions file." << std::endl;
    std::cout << "\t-t: (Optional) Run CFL check after simulation." << std::endl;
}

static bool get_args(int argc, char *argv[], bool *write_test, bool *run_cfl, std::string *cfg_path, std::string *initial_conds_path) {
    int ch = 0;
    while ((ch = getopt(argc, argv, "wtc:s:")) != -1) {
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
            case 't':
                *run_cfl = true;
                break;
            default:
                return false;
        }
    }
    return true;
}

template<typename T>
requires Numeric<T>
void run_simulation_internal(const SimulationConfig &config, const std::string &initial_conds_path, bool run_cfl) {
    auto initial_conditions = read_initial_conditions<T>(initial_conds_path);
    auto flux = match_flux<T>(config.flux);
    // For now, only difference solvers.
    auto solver = match_difference<T>(config.solver);
    auto solution = solver->solve(initial_conditions, config.discretization_size, config.num_timesteps, config.delta_t, config.delta_x, flux);
    solution.print_system();

    if (run_cfl) {
        solver->cfl_check_mesh(solution, flux, config.delta_t, config.delta_x);
    }

    delete solver;
    delete flux;
}

void run_simulation(const std::string &cfg_path, const std::string &initial_conds_path, bool run_cfl) {
    // Read config
    auto config = read_config(cfg_path);

    // Note: duplicated logic between these branches because we can only initialize once we know the templated type.
    if (config.domain == "real") {
        run_simulation_internal<Real>(config, initial_conds_path, run_cfl);
    } else if (config.domain == "interval") {
        run_simulation_internal<Winterval>(config, initial_conds_path, run_cfl);
    } else if (config.domain == "affine") {
        run_simulation_internal<AffineForm>(config, initial_conds_path, run_cfl);
    } else if (config.domain == "mixed") {
        run_simulation_internal<MixedForm>(config, initial_conds_path, run_cfl);
    } else {
        std::cerr << "Invalid domain!" << std::endl;
        exit(EXIT_FAILURE);
    }
}
