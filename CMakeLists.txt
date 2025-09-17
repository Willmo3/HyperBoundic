cmake_minimum_required(VERSION 4.0)
project(PDEapprox)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(Winterval)

add_executable(PDEapprox src/main.cpp
        src/domains/Real.cpp
        src/domains/Real.hpp
        src/domains/Numeric.hpp
        src/solvers/LaxFriedrichsSolver.hpp
        src/discretizations/PdeDiscretization.hpp
)