include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# Windows compatibility
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../out/tests)

# Scheme tests
add_executable(test_friedrichs difference/test_friedrichs.cpp)
target_link_libraries(test_friedrichs GTest::gtest_main)
add_executable(test_leapfrog difference/test_leapfrog.cpp)
target_link_libraries(test_leapfrog GTest::gtest_main)
add_executable(test_serialization difference/test_serialization.cpp)
target_link_libraries(test_serialization GTest::gtest_main)

# Flux tests
add_executable(test_flux difference/test_flux.cpp)
target_link_libraries(test_flux GTest::gtest_main)

include(GoogleTest)
target_link_libraries(test_friedrichs friedrichs)
target_link_libraries(test_leapfrog leapfrog)
target_link_libraries(test_serialization friedrichs)
target_link_libraries(test_flux friedrichs)

# Visualization executables
add_executable(visualize_leapfrog viz/visualize_leapfrog.cpp)
target_link_libraries(visualize_leapfrog leapfrog matplot)

add_executable(visualize_friedrichs viz/visualize_friedrichs.cpp)
target_link_libraries(visualize_friedrichs friedrichs matplot)
