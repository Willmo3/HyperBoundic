cmake_minimum_required(VERSION 3.5)
project(PDEnclose)

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 20)
include_directories(lib)
include_directories(src)

add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(test)
