//
// Created by will on 11/3/25.
//

#ifndef PDENCLOSE_READ_DISCRETIZATIONS_H
#define PDENCLOSE_READ_DISCRETIZATIONS_H
#include <fstream>
#include <string>

#include "cereal/archives/json.hpp"
#include "cereal/types/vector.hpp"
#include "domains/Numeric.hpp"

/**
 *
 * @tparam T Numeric type of initial conditions.
 * @param file_name Name of file to read initial conditions from.
 * @return A vector of the initial conditions, read from file_name.
 */
template<typename T>
requires Numeric<T>
std::vector<T> read_initial_conditions(const std::string& file_name) {
    std::ifstream f;
    f.open(file_name);

    auto initial_conds = std::vector<T>();
    {
        cereal::JSONInputArchive archive(f);
        archive(initial_conds);
    }

    f.close();
    return initial_conds;
}

/**
 *
 * @tparam T Numeric type of initial conditions
 * @param file_name name of file to write conditions to.
 * @param initial_conds Initial conditions to write
 */
template<typename T>
requires Numeric<T>
void write_initial_conditions(const std::string& file_name, const std::vector<T> &initial_conds) {
    std::ofstream f;
    f.open(file_name);
    {
        cereal::JSONOutputArchive archive(f);
        archive(initial_conds);
    }
    f.close();
}

#endif //PDENCLOSE_READ_DISCRETIZATIONS_H
