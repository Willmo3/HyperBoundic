//
// Created by will on 11/7/25.
//

#include "generate_initial_conditions.hpp"

#include "Caffeine/AffineForm.hpp"
#include "domains/Real.hpp"
#include "DualDomain/MixedForm.hpp"
/*
 * Conversion function.
 * Convert original conditions to other types with tolerance 0.1.
 */
std::vector<Winterval> convert_conds_to_interval(const std::vector<Real> &real_conds, double epsilon) {
    epsilon = std::abs(epsilon);
    std::vector<Winterval> interval_conds = std::vector<Winterval>(real_conds.size());

    for (int i = 0; i < real_conds.size(); i++) {
        interval_conds[i] = Winterval(real_conds[i].value() - epsilon, real_conds[i].value() + epsilon);
    }
    return interval_conds;
}
std::vector<AffineForm> convert_conds_to_affine(const std::vector<Real> &real_conds, double epsilon) {
    epsilon = std::abs(epsilon);
    std::vector<AffineForm> affine_conds = std::vector<AffineForm>(real_conds.size());

    for (int i = 0; i < real_conds.size(); i++) {
        affine_conds[i] = AffineForm(Winterval(real_conds[i].value() - epsilon, real_conds[i].value() + epsilon));
    }
    return affine_conds;
}
std::vector<MixedForm> convert_conds_to_mixed(const std::vector<Real> &real_conds, double epsilon) {
    epsilon = std::abs(epsilon);
    std::vector<MixedForm> mixed_conds = std::vector<MixedForm>(real_conds.size());

    for (int i = 0; i < real_conds.size(); i++) {
        mixed_conds[i] = MixedForm(Winterval(real_conds[i].value() - epsilon, real_conds[i].value() + epsilon));
    }
    return mixed_conds;
}

/*
 * Write initial conditions for different fluxes.
 */
void write_burgers_initial_conds(const std::string &root) {
    auto size = 20;
    auto tolerance = 0.05;

    // Prepare base initial conditons.
    std::vector<Real> base_conds = std::vector<Real>(size);
    for (int i = 0; i < base_conds.size(); i++) {
        auto x = i * 2;
        base_conds[i] = x < 15.01 ? -0.015 * x * (x - 15) : 0;
    }
    write_initial_conditions<Real>(root + "/burgers_real_conds.json", base_conds);

    auto interval_conds = convert_conds_to_interval(base_conds, tolerance);
    write_initial_conditions<Winterval>(root + "/burgers_interval_conds.json", interval_conds);

    auto affine_conds = convert_conds_to_affine(base_conds, tolerance);
    write_initial_conditions<AffineForm>(root + "/burgers_affine_conds.json", affine_conds);

    auto mixed_conds = convert_conds_to_mixed(base_conds, tolerance);
    write_initial_conditions<MixedForm>(root + "/burgers_mixed_conds.json", mixed_conds);
}
void write_lwr_initial_conds(const std::string &root) {
    auto size = 20;
    auto tolerance = 0.07;

    // Prepare base initial conditons.
    std::vector<Real> base_conds = std::vector<Real>(size);
    for (int i = 0; i < base_conds.size(); i++) {
        auto x = i * 2;
        base_conds[i] = x < 15.01 ? -0.015 * x * (x - 15) : 0;
    }
    write_initial_conditions<Real>(root + "/lwr_real_conds.json", base_conds);

    auto interval_conds = convert_conds_to_interval(base_conds, tolerance);
    write_initial_conditions<Winterval>(root + "/lwr_interval_conds.json", interval_conds);

    auto affine_conds = convert_conds_to_affine(base_conds, tolerance);
    write_initial_conditions<AffineForm>(root + "/lwr_affine_conds.json", affine_conds);

    auto mixed_conds = convert_conds_to_mixed(base_conds, tolerance);
    write_initial_conditions<MixedForm>(root + "/lwr_mixed_conds.json", mixed_conds);
}
void write_buckley_leverett_initial_conds(const std::string &root) {
    auto size = 20;
    auto tolerance = 0.05;

    // Prepare base initial conditons.
    std::vector<Real> base_conds = std::vector<Real>(size);
    for (int i = 0; i < base_conds.size(); i++) {
        auto x = i * 2;
        base_conds[i] = x < 15.01 ? -0.015 * x * (x - 15) : 0;
    }
    write_initial_conditions<Real>(root + "/buckley_leverett_real_conds.json", base_conds);

    auto interval_conds = convert_conds_to_interval(base_conds, tolerance);
    write_initial_conditions<Winterval>(root + "/buckley_leverett_interval_conds.json", interval_conds);

    auto affine_conds = convert_conds_to_affine(base_conds, tolerance);
    write_initial_conditions<AffineForm>(root + "/buckley_leverett_affine_conds.json", affine_conds);

    auto mixed_conds = convert_conds_to_mixed(base_conds, tolerance);
    write_initial_conditions<MixedForm>(root + "/buckley_leverett_mixed_conds.json", mixed_conds);
}
void write_cubic_initial_conds(const std::string &root) {
    auto size = 20;
    auto tolerance = 0.1;

    // Prepare base initial conditons.
    std::vector<Real> base_conds = std::vector<Real>(size);
    for (int i = 0; i < base_conds.size(); i++) {
        auto x = i * 2;
        base_conds[i] = x < 15.01 ? -0.015 * x * (x - 15) : 0;
    }
    write_initial_conditions<Real>(root + "/cubic_real_conds.json", base_conds);

    auto interval_conds = convert_conds_to_interval(base_conds, tolerance);
    write_initial_conditions<Winterval>(root + "/cubic_interval_conds.json", interval_conds);

    auto affine_conds = convert_conds_to_affine(base_conds, tolerance);
    write_initial_conditions<AffineForm>(root + "/cubic_affine_conds.json", affine_conds);

    auto mixed_conds = convert_conds_to_mixed(base_conds, tolerance);
    write_initial_conditions<MixedForm>(root + "/cubic_mixed_conds.json", mixed_conds);
}

// write initial conditions for the domains
void generate_initial_conds(const std::string& root) {
    // TODO: ask: is there a notion of meaningful initial conditions for certain tasks?
    auto size = 20;
    // Prepare base initial conditons.
    std::vector<Real> base_conds = std::vector<Real>(size);
    for (int i = 0; i < base_conds.size(); i++) {
        base_conds[i] = i * 0.2 + 1;
    }
    write_lwr_initial_conds(root);
    write_buckley_leverett_initial_conds(root);
    write_cubic_initial_conds(root);
    write_burgers_initial_conds(root);
}