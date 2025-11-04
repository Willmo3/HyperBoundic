add_library(domains
        domains/Real.cpp
        domains/Real.hpp
        domains/Numeric.hpp
)
target_link_libraries(domains winterval caffeine dualdomain)

add_library(fluxes
        solvers/flux/FluxFunction.hpp
        solvers/flux/CubicFlux.hpp
        solvers/flux/BurgersFlux.hpp
        solvers/flux/LwrFlux.hpp
        solvers/flux/BuckleyLeverettFlux.hpp
        solvers/flux/match_flux.hpp
        domains/Numeric.hpp
)
set_target_properties(fluxes PROPERTIES LINKER_LANGUAGE CXX)

add_library(discretizations
        meshes/RectangularMesh.hpp
        meshes/CflCheck.hpp
        domains/Numeric.hpp
)
set_target_properties(discretizations PROPERTIES LINKER_LANGUAGE CXX)

# Simulations
add_library(friedrichs solvers/difference/LaxFriedrichsSolver.hpp)
target_link_libraries(friedrichs domains fluxes discretizations)

add_library(leapfrog solvers/difference/LeapfrogSolver.hpp)
target_link_libraries(leapfrog domains fluxes discretizations)

add_library(local_friedrichs solvers/volume/LocalLaxFriedrichsSolver.hpp)
target_link_libraries(local_friedrichs domains fluxes discretizations)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../out)

add_executable(PDEapprox
        exe/main.cpp
        exe/initial_conditions.hpp
        exe/SimulationConfig.hpp
        visualization/MeshVisualizer.hpp)
target_link_libraries(PDEapprox friedrichs leapfrog local_friedrichs)
target_link_libraries(PDEapprox matplot)