cmake_minimum_required(VERSION 4.0)
project(PDEapprox)

set(CMAKE_CXX_STANDARD 20)

add_executable(PDEapprox main.cpp
        src/PdeDiscretization.cpp
        src/PdeDiscretization.h
        src/LaxFriedrichsSolver.cpp
        src/LaxFriedrichsSolver.h)
