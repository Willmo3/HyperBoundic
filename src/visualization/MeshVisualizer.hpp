//
// Created by will on 9/20/25.
//

#ifndef PDEAPPROX_DISCRETIZATIONVISUALIZER_H
#define PDEAPPROX_DISCRETIZATIONVISUALIZER_H
#include "meshes/RectangularMesh.hpp"
#include "matplot/freestanding/plot.h"
#include "matplot/util/common.h"

/*
 * Note: we use individual functions, as opposed to a shared class,
 * to allow the transform lambda to capture the system.
 */
inline void show_real_surface(RectangularMesh<Real> *system) {
    using namespace matplot;
    auto [X, Y] = meshgrid(iota(0, 1, 3));

    auto Z = transform(X, Y, [system](double x_coord, double y_coord) {
        auto x_int = static_cast<uint32_t>(x_coord);
        auto y_int = static_cast<uint32_t>(y_coord);
        return system->get(x_int, y_int).value();
    });

    auto data = surfc(X, Y, Z);
    xlabel("timestep");
    ylabel("space");
    zlabel("wave height");
    show();
}

#endif //PDEAPPROX_DISCRETIZATIONVISUALIZER_H