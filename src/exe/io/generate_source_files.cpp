//
// Created by will on 11/4/25.
//

#include <filesystem>
#include <iostream>

#include "initial_conditions.hpp"
#include "SimulationConfig.hpp"
#include "Caffeine/AffineForm.hpp"
#include "domains/Real.hpp"
#include "DualDomain/MixedForm.hpp"
#include "Winterval/Winterval.hpp"

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

// write initial conditions for the domains
void write_conditions(const std::string& root) {
    auto size = 20;
    std::vector<Real> real_conds = std::vector<Real>(size);

    for (int i = 0; i < real_conds.size(); i++) {
        real_conds[i] = i * 0.2 + 1;
    }
    write_initial_conditions<Real>(root + "/real_conds.json", real_conds);

    std::vector<Winterval> interval_conds = std::vector<Winterval>(size);
    for (int i = 0; i < interval_conds.size(); i++) {
        interval_conds[i] = Winterval(i * 0.2 + 0.9, i * 0.2 + 1.1);
    }
    write_initial_conditions<Winterval>(root + "/interval_conds.json", interval_conds);

    std::vector<AffineForm> affine_conds = std::vector<AffineForm>(size);
    for (int i = 0; i < affine_conds.size(); i++) {
        affine_conds[i] = AffineForm(Winterval(i * 0.2 + 0.9, i * 0.2 + 1.1));
    }
    write_initial_conditions<AffineForm>(root + "/affine_conds.json", affine_conds);

    std::vector<MixedForm> mixed_conds = std::vector<MixedForm>(size);
    for (int i = 0; i < mixed_conds.size(); i++) {
        mixed_conds[i] = MixedForm(Winterval(i * 0.2 + 0.9, i * 0.2 + 1.1));
    }
    write_initial_conditions<MixedForm>(root + "/mixed_conds.json", mixed_conds);
}

// Write out all the config files for the various combinations.
void write_config_files(const std::string &root) {
    std::string domains[] = {"real", "interval", "affine", "mixed"};
    std::string fluxes[] = {"cubic", "burgers", "lwr", "buckley_leverett"};
    std::string solvers[] = {"lax_friedrichs", "leapfrog" };

    for (auto domain: domains) {
        for (auto flux: fluxes) {
            for (auto solver: solvers) {
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

    write_conditions(path);
    write_config_files(path);
}