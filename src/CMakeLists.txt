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

# Serialization
add_library(io io/save_initial_conditions.hpp)
target_link_libraries(io domains)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../out)

add_executable(PDEapprox main.cpp visualization/MeshVisualizer.hpp)
target_link_libraries(PDEapprox friedrichs leapfrog local_friedrichs io)
target_link_libraries(PDEapprox matplot)