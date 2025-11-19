include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# Windows compatibility
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../out/tests)

# Difference tests
add_executable(test_friedrichs difference/test_friedrichs.cpp)
target_link_libraries(test_friedrichs GTest::gtest_main)
add_executable(test_leapfrog difference/test_leapfrog.cpp)
target_link_libraries(test_leapfrog GTest::gtest_main)
add_executable(test_serialization difference/test_serialization.cpp)
target_link_libraries(test_serialization GTest::gtest_main)

# Volume tests
add_executable(test_local_lax_friedrichs volume/test_local_friedrichs.cpp)
target_link_libraries(test_local_lax_friedrichs GTest::gtest_main)

# Flux tests
add_executable(test_flux difference/test_flux.cpp)
target_link_libraries(test_flux GTest::gtest_main)

include(GoogleTest)
target_link_libraries(test_friedrichs difference_solvers)
target_link_libraries(test_leapfrog difference_solvers)
target_link_libraries(test_serialization difference_solvers)
# we only test flux functions w/ difference meshes bc it makes no difference on underlying math
target_link_libraries(test_flux difference_solvers)
target_link_libraries(test_local_lax_friedrichs volume_solvers)

# Visualization executables
add_executable(visualize_leapfrog viz/visualize_leapfrog.cpp)
target_link_libraries(visualize_leapfrog difference_solvers matplot)

add_executable(visualize_friedrichs viz/visualize_friedrichs.cpp)
target_link_libraries(visualize_friedrichs difference_solvers matplot)
