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
    std::vector<Real> real_conds = std::vector<Real>(4);
    real_conds[0] = 1;
    real_conds[1] = 2;
    real_conds[2] = 3;
    real_conds[3] = 4;
    write_initial_conditions<Real>(root + "/real_conds.json", real_conds);

    std::vector<Winterval> interval_conds = std::vector<Winterval>(4);
    interval_conds[0] = Winterval(0.9, 1.1);
    interval_conds[1] = Winterval(1.9, 2.1);
    interval_conds[2] = Winterval(2.9, 3.1);
    interval_conds[3] = Winterval(3.9, 4.1);
    write_initial_conditions<Winterval>(root + "/interval_conds.json", interval_conds);

    std::vector<AffineForm> affine_conds = std::vector<AffineForm>(4);
    affine_conds[0] = AffineForm(Winterval(0.9, 1.1));
    affine_conds[1] = AffineForm(Winterval(1.9, 2.1));
    affine_conds[2] = AffineForm(Winterval(2.9, 3.1));
    affine_conds[3] = AffineForm(Winterval(3.9, 4.1));
    write_initial_conditions<AffineForm>(root + "/affine_conds.json", affine_conds);

    std::vector<MixedForm> mixed_conds = std::vector<MixedForm>(4);
    mixed_conds[0] = MixedForm(Winterval(0.9, 1.1));
    mixed_conds[1] = MixedForm(Winterval(1.9, 2.1));
    mixed_conds[2] = MixedForm(Winterval(2.9, 3.1));
    mixed_conds[3] = MixedForm(Winterval(3.9, 4.1));
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
                auto cfg = SimulationConfig(domain, flux, solver, 4, 4, 1, 0.01);
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