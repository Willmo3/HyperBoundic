add_library(lf
        domains/Real.cpp
        domains/Real.hpp
        domains/Numeric.hpp
        solvers/LaxFriedrichsSolver.hpp
        discretizations/PdeDiscretization.hpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../out)

add_executable(PDEapprox main.cpp)
target_link_libraries(PDEapprox lf)
target_link_libraries(PDEapprox winterval)

