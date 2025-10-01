cmake_minimum_required(VERSION 3.5)
project(PDEapprox)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(lib/Waffine)
add_subdirectory(lib/Winterval)
add_subdirectory(lib/matplotplusplus)
add_subdirectory(src)
add_subdirectory(test)
