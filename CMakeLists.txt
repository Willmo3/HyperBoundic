cmake_minimum_required(VERSION 3.5)
project(PDEnclose)

set(CMAKE_CXX_STANDARD 20)
include_directories(lib)
include_directories(src)

find_package(OpenMP REQUIRED)

add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(test)
