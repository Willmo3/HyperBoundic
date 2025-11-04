add_library(domains
        domains/Real.cpp
        domains/Real.hpp
        domains/Numeric.hpp
)
target_link_libraries(domains winterval caffeine dualdomain)

add_library(fluxes
        flux/FluxFunction.hpp
        flux/CubicFlux.hpp
        flux/BurgersFlux.hpp
        flux/LwrFlux.hpp
        flux/BuckleyLeverettFlux.hpp
        flux/match_flux.hpp
        domains/Numeric.hpp
)
set_target_properties(fluxes PROPERTIES LINKER_LANGUAGE CXX)

add_library(discretizations
        meshes/RectangularMesh.hpp
        meshes/CflCheck.hpp
        domains/Numeric.hpp
)
set_target_properties(discretizations PROPERTIES LINKER_LANGUAGE CXX)

add_library(solvers
        solvers/difference/LaxFriedrichsSolver.hpp
        solvers/difference/LeapfrogSolver.hpp
        solvers/volume/LocalLaxFriedrichsSolver.hpp
        solvers/PDESolver.hpp
)
target_link_libraries(solvers domains fluxes discretizations)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../out)

add_executable(PDEapprox
        exe/main.cpp
        exe/initial_conditions.hpp
        exe/SimulationConfig.hpp
        visualization/MeshVisualizer.hpp)
target_link_libraries(PDEapprox solvers matplot)