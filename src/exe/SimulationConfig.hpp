//
// Created by will on 11/3/25.
//

#ifndef PDENCLOSE_SIMULATIONCONFIG_H
#define PDENCLOSE_SIMULATIONCONFIG_H
#include <fstream>
#include <string>
#include <utility>

#include "cereal/details/helpers.hpp"

#include "cereal/archives/json.hpp"

/**
 * Config options for PDE simulation.
 */
struct SimulationConfig {
    /*
     * Constructors
     */

    /**
     * @param flux Name of flux function being serialized. Options: cubic, burgers, lwr, buckley_leverett
     * @param domain Name of abstract domain serialized over. Options: real, interval, affine, mixed
     * @param discretization_size Size of the discretization being serialized.
     * @param num_timesteps Number of timesteps to run simulation for.
     * @param delta_x spatial step
     * @param delta_t timestep
     */
    SimulationConfig(
        std::string flux, std::string domain, uint32_t discretization_size, uint32_t num_timesteps, double delta_x, double delta_t):
          flux(std::move(flux)),
          domain(std::move(domain)),
          discretization_size(discretization_size),
          num_timesteps(num_timesteps),
          delta_t(delta_t), delta_x(delta_x) {}

    /**
     * Empty default constructor allows serialization.
     */
    SimulationConfig() = default;

    template<class Archive>
    void serialize(Archive & archive) {
        archive(cereal::make_nvp("flux", flux),
                cereal::make_nvp("domain", domain),
                cereal::make_nvp("discretization_size", discretization_size),
                cereal::make_nvp("timesteps", num_timesteps),
                cereal::make_nvp("delta_x", delta_x),
                cereal::make_nvp("delta_t", delta_t));
    }

    /*
     * Fields
     */
    std::string flux;
    std::string domain;
    uint32_t discretization_size;
    uint32_t num_timesteps;
    double delta_t;
    double delta_x;
};

/*
 * io
 */

/**
 * @brief Read a configuration tuple from an input file.
 * This function is static to avoid lifetime issues with Cereal.
 * @param file_name Name of file to read config in from.
 * @return A config tuple for a simulation.
 */
SimulationConfig read_config(const std::string &file_name) {
    std::ifstream f;
    f.open(file_name);

    auto config = SimulationConfig();
    {
        cereal::JSONInputArchive archive(f);
        archive(config);
    }

    f.close();
    return config;
}

/**
 * Write a configuration tuple to an output file.
 * This function is static to avoid lifetime issues with Cereal.
 *
 * @param file_name Name of file to write configuration to.
 * @param config Configuration to write.
 */
void write_config(const std::string &file_name, const SimulationConfig &config) {
    std::ofstream f;
    f.open(file_name);
    {
        cereal::JSONOutputArchive archive(f);
        archive(config);
    }
    f.close();
}


#endif //PDENCLOSE_SIMULATIONCONFIG_H
