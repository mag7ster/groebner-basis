cmake_minimum_required(VERSION 3.27.4)

project(groebner_basis)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -Wall")

set(BOOST_INCLUDE_LIBRARIES rational)
add_subdirectory(deps/boost EXCLUDE_FROM_ALL)

set(SOURCE_EXE main.cpp)
set(GROEBNER_LIB groebner_basis.cpp)
set(HELP_LIB tostringhelper.cpp)


add_library(groebner_basis STATIC ${GROEBNER_LIB})
add_library(tostringhelper STATIC ${HELP_LIB})
add_executable(main ${SOURCE_EXE})

target_link_libraries(main groebner_basis Boost::rational)
target_link_libraries(groebner_basis tostringhelper)
target_link_libraries(tostringhelper Boost::rational)

