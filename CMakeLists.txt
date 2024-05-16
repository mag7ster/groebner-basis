cmake_minimum_required(VERSION 3.27.4)

project(groebner_basis)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(TESTS_EXE src/tests.cpp)
set(BENCH_EXE src/bench.cpp)

set(Boost_USE_STATIC_LIBS        ON)  # only find static libs
set(Boost_USE_DEBUG_LIBS        OFF)  # ignore debug libs and
set(Boost_USE_RELEASE_LIBS       ON)  # only find release libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     ON)  # link Boost Static libraries

find_package(Boost)
include_directories( ${Boost_INCLUDE_DIRS} )

enable_testing()
find_package(GTest CONFIG REQUIRED)
include(GoogleTest)
add_executable(tests ${TESTS_EXE})
target_link_libraries(tests PRIVATE GTest::gtest_main GTest::gtest)
add_test(AllTestsInMain tests)
gtest_discover_tests(tests)

find_package(benchmark CONFIG REQUIRED)
add_executable(bench ${BENCH_EXE})
target_link_libraries(bench PRIVATE benchmark::benchmark benchmark::benchmark_main)
