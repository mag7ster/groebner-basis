cmake_minimum_required(VERSION 3.27.4)

project(groebner_basis)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(SOURCE_EXE src/main.cpp)
# set(GROEBNER_LIB src/groebner_basis.cpp)

set(Boost_USE_STATIC_LIBS        ON)  # only find static libs
set(Boost_USE_DEBUG_LIBS        OFF)  # ignore debug libs and
set(Boost_USE_RELEASE_LIBS       ON)  # only find release libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     ON)  # link Boost Static libraries

find_package(Boost)
include_directories( ${Boost_INCLUDE_DIRS} )

# add_library(groebner_basis STATIC ${GROEBNER_LIB})
add_executable(main ${SOURCE_EXE})
