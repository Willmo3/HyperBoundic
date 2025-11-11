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
        domains/Numeric.hpp
)
set_target_properties(fluxes PROPERTIES LINKER_LANGUAGE CXX)

add_library(discretizations
        meshes/RectangularMesh.hpp
        meshes/CflCheck.hpp
        domains/Numeric.hpp
)
set_target_properties(discretizations PROPERTIES LINKER_LANGUAGE CXX)

add_library(difference_solvers
        solvers/difference/LaxFriedrichsSolver.hpp
        solvers/difference/LeapfrogSolver.hpp
        solvers/difference/DifferenceSolver.hpp
)
target_link_libraries(difference_solvers domains fluxes discretizations)
add_library(volume_solvers
        solvers/volume/LocalLaxFriedrichsSolver.hpp
)
target_link_libraries(volume_solvers domains fluxes discretizations)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../out)

add_executable(PDEapprox
        exe/main.cpp
        exe/experiment/SimulationConfig.hpp
        exe/experiment/generators/generate_source_files.cpp
        exe/experiment/generators/generate_source_files.h
        exe/experiment/generators/generate_config_files.cpp
        exe/experiment/generators/generate_config_files.h
        exe/experiment/generators/generate_initial_conditions.cpp
        exe/experiment/generators/generate_initial_conditions.hpp
        exe/args/match_names.hpp
        visualization/MeshVisualizer.hpp)
target_link_libraries(PDEapprox difference_solvers volume_solvers matplot)