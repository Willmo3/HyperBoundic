add_library(domains
        domains/Real.cpp
        domains/Real.hpp
        domains/Numeric.hpp
)
target_link_libraries(domains winterval waffine wixed)

add_library(lf
        solvers/LaxFriedrichsSolver.hpp
        solvers/LeapfrogSolver.hpp
        solvers/flux/FluxFunction.hpp
        solvers/flux/CubicFlux.hpp
        solvers/cfl/CflCheck.hpp
        discretizations/PdeDiscretization.hpp
)
target_link_libraries(lf domains)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../out)

add_executable(PDEapprox main.cpp visualization/DiscretizationVisualizers.hpp)
target_link_libraries(PDEapprox lf)
target_link_libraries(PDEapprox matplot)

