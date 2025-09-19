add_library(lf
        domains/Real.cpp
        domains/Real.hpp
        domains/Numeric.hpp
        solvers/LaxFriedrichsSolver.hpp
        solvers/flux/FluxFunction.hpp
        solvers/flux/CubicFlux.hpp
        solvers/cfl/cfl_check.hpp
        discretizations/PdeDiscretization.hpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../out)

add_executable(PDEapprox main.cpp)
target_link_libraries(PDEapprox lf)
target_link_libraries(PDEapprox winterval)
target_link_libraries(PDEapprox matplot)

