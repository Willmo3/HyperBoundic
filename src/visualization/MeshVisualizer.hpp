//
// Created by will on 9/20/25.
//

#ifndef PDENCLOSE_DISCRETIZATIONVISUALIZER_H
#define PDENCLOSE_DISCRETIZATIONVISUALIZER_H
#include "meshes/RectangularMesh.hpp"
#include "matplot/freestanding/plot.h"
#include "matplot/util/common.h"

/*
 * Note: we use individual functions, as opposed to a shared class,
 * to allow the transform lambda to capture the system.
 */
inline void prepare_matplot(RectangularMesh<Real> *system) {
    using namespace matplot;
    auto [X, Y] =
        meshgrid(iota(0, 1, std::min(system->discretization_size() - 1, system->num_timesteps() - 1)));

    auto Z = transform(X, Y, [system](double x_coord, double y_coord) {
        auto x_int = static_cast<uint32_t>(x_coord);
        auto y_int = static_cast<uint32_t>(y_coord);
        return system->get(x_int, y_int).value();
    });

    auto data = ribbon(X, Y, Z);
    axis({0, static_cast<double>(system->num_timesteps()), 0, static_cast<double>(system->discretization_size())});
    xlabel("timestep");
    ylabel("space");
    zlabel("wave height");
    colorbar();
}

inline void show_real_surface(RectangularMesh<Real> *system) {
    using namespace matplot;
    prepare_matplot(system);
    show();
}

inline void save_real_surface(RectangularMesh<Real> *system, const std::string &filename) {
    using namespace matplot;
    prepare_matplot(system);
    save(filename);
}

#endif //PDENCLOSE_DISCRETIZATIONVISUALIZER_H