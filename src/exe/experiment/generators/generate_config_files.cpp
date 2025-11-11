//
// Created by will on 11/7/25.
//

#include "generate_config_files.h"

#include "exe/experiment/SimulationConfig.hpp"

void write_single_cfg(const std::string &domain_name, const std::string &solver_name, const std::string &flux_name, double timestep, const std::string &root) {
    // for now, these are constant between all simulations
    auto discretization_size = 20;
    auto num_timesteps = 25;
    auto delta_x = 2;

    auto cfg = SimulationConfig(domain_name, flux_name, solver_name, discretization_size, num_timesteps, delta_x, timestep);
    auto file_name = root + "/" + domain_name + "_" + flux_name + "_" + solver_name + "_config.json";
    write_config(file_name, cfg);
}

void generate_config_files(const std::string &root) {
    std::string domains[] = {"real", "interval", "affine", "mixed"};
    std::string fluxes[] = {"cubic", "burgers", "lwr", "buckley_leverett"};
    std::string solvers[] = {"lax_friedrichs", "leapfrog" };

    for (const auto& domain: domains) {
        for (const auto& solver: fluxes) {
            auto burgers_timestep = 1;
            auto lwr_timestep = 25;
            // Cubic is highly nonlinear. For now, we're going to use a smaller timestep than Phocus.
            // Note that Phocus used different timesteps for the different solvers here.
            auto cubic_timestep = 0.01;
            auto buckley_leverett_timestep = 0.35;

            // Generate timestepped config files for each flux function.
            write_single_cfg(domain, solver, "burgers", burgers_timestep, root);
            write_single_cfg(domain, solver, "lwr", lwr_timestep, root);
            write_single_cfg(domain, solver, "cubic", cubic_timestep, root);
            write_single_cfg(domain, solver, "buckley_leverett", buckley_leverett_timestep, root);
        }
    }
}
