//
// Created by will on 11/4/25.
//

#include <filesystem>
#include <iostream>

#include "generate_initial_conditions.hpp"
#include "../SimulationConfig.hpp"

// Create a simulations directory in the current working directory if it does not exist.
void create_dir() {
    namespace fs = std::filesystem;
    fs::path dir = "simulations";

    std::error_code ec;
    if (fs::exists(dir, ec)) {
        if (!fs::is_directory(dir, ec)) {
            std::cerr << "Path " << dir << " exists but is not a directory\n";
        }
        // directory already exists -> nothing to do
        return;
    }

    if (!fs::create_directories(dir, ec)) {
        std::cerr << "Failed to create directory " << dir << ": " << ec.message() << '\n';
    }
}

// Write out all the config files for the various combinations.
void generate_config_files(const std::string &root) {
    std::string domains[] = {"real", "interval", "affine", "mixed"};
    std::string fluxes[] = {"cubic", "burgers", "lwr", "buckley_leverett"};
    std::string solvers[] = {"lax_friedrichs", "leapfrog" };
    // TODO: configure different timesteps for each

    for (const auto& domain: domains) {
        for (const auto& flux: fluxes) {
            for (const auto& solver: solvers) {
                auto cfg = SimulationConfig(domain, flux, solver, 20, 25, 1, 0.01);
                auto file_name = root + "/" + domain + "_" + flux + "_" + solver + "_config.json";
                write_config(file_name, cfg);
            }
        }
    }
}

void generate_source_files() {
    create_dir();
    auto path = std::filesystem::current_path().string() + "/simulations";
    std::cout << "Generated source files will be placed in: " << path << std::endl;

    generate_initial_conds(path);
    generate_config_files(path);
}