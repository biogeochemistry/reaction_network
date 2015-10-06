########################################################################
# CMake build script for Google Test.
#
# To run the tests for Google Test itself on Linux, use 'make test' or
# ctest.  You can select which tests to run using 'ctest -R regex'.
# For more options, run 'ctest --help'.

# BUILD_SHARED_LIBS is a standard CMake variable, but we declare it here to
# make it prominent in the GUI.
option(BUILD_SHARED_LIBS "Build shared libraries (DLLs)." ON)

# When other libraries are using a shared version of runtime libraries,
# Google Test also has to use one.
option(  gtest_force_shared_crt
  "Use shared (DLL) run-time lib even when Google Test is built as static lib."
  ON)

option(test "Build all of gtest's own tests." ON)


# Defines pre_project_set_up_hermetic_build() and set_up_hermetic_build().
include(cmake/hermetic_build.cmake OPTIONAL)

if (COMMAND pre_project_set_up_hermetic_build)
  pre_project_set_up_hermetic_build()
endif()

########################################################################
#
# Project-wide settings

# Name of the project.
#
# CMake files in this project can refer to the root source directory
# as ${gtest_SOURCE_DIR} and to the root binary directory as
# ${gtest_BINARY_DIR}.
# Language "C" is required for find_package(Threads).
project(reaction_network CXX)
cmake_minimum_required(VERSION 2.6.2)

if (COMMAND set_up_hermetic_build)
  set_up_hermetic_build()
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)

find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)

# compiler options:
# set(CMAKE_CXX_COMPILER "clang++") #compiler clang++
# set(CMAKE_CXX_COMPILER "g++") #compiler g++
# set(CMAKE_CXX_COMPILER "/usr/local/Cellar/gcc49/4.9.0/bin/g++-4.9") #compiler g++
set (CMAKE_CXX_LINK_FLAGS "-fopenmp -lstdc++")
add_definitions("-std=c++11 -O3 -ffast-math -ffinite-math-only -fopenmp -fPIC -Wall -Wextra -pthread -Wno-missing-field-initializers  -fdiagnostics-color=always")
# add_definitions("-std=c++11 -g -O0 -fopenmp -fPIC -Wall -Wextra -Wno-missing-field-initializers  -fdiagnostics-color=always -pthread")

# CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# folders to import
set(reaction_network ${CMAKE_SOURCE_DIR})
set(GTEST_DIR "${reaction_network}/libs/gtest/")
set(EIGEN_DIR "${reaction_network}/libs/eigen/")
set(GNUPLOT_DIR "${reaction_network}/libs/gnuplot/")
set(GNUPLOT2_DIR "${reaction_network}/libs/gnuplot-iostream/")
set(INCLUDE_DIR "${reaction_network}/include/")


# Define helper functions and macros used by Google Test.
include(${GTEST_DIR}/cmake/internal_utils.cmake)

# # config_compiler_and_linker()  # Defined in internal_utils.cmake.

# # Where Google Test's .h files can be found.
include_directories(
  ${GTEST_DIR}/include
  ${GTEST_DIR}
  ${EIGEN_DIR}
  ${GNUPLOT_DIR}
  ${GNUPLOT2_DIR}
  ${CMAKE_SOURCE_DIR}
  ${INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS})

# # Where Google Test's libraries can be found.
link_directories(${GTEST_DIR}/src)

set(HEADER_FILES  
                  include/OdeSolver.hpp 
                  include/Runge-Kutta.hpp
                  include/BoundaryConditions.hpp 
                  include/SecondOrderOde.hpp 
                  include/FiniteDifferenceGrid.hpp 
                  include/BvpOde.hpp 
                  include/Node.hpp 
                  include/LinearSolver.hpp
                  include/BvpPde.hpp
                  include/CoupledPde.hpp)
set(SOURCES_FILES 
                  src/OdeSolver.cpp
                  src/Runge-Kutta.cpp 
                  src/BoundaryConditions.cpp  
                  src/FiniteDifferenceGrid.cpp 
                  src/LinearSolver.cpp 
                  src/BvpOde.cpp 
                  src/SecondOrderOde.cpp
                  src/BvpPde.cpp
                  src/CoupledPde.cpp)

# ########################################################################
# #
# # Defines the gtest & gtest_main libraries.  User tests should link
# # with one of them.

# # Google Test libraries.  We build them using more strict warnings than what
# # are used for other targets, to ensure that gtest can be compiled by a user
# # aggressive about warnings.

cxx_library(gtest "${cxx_strict}" ${GTEST_DIR}/src/gtest-all.cc)
cxx_library(gtest_main "${cxx_strict}" ${GTEST_DIR}/src/gtest_main.cc)


# # cxx_executable(name dir lib srcs...)
# #
# # creates a named target that depends on the given lib and is built
# # from the given source files.  dir/name.cc is implicitly included in
# # the source file list.

ADD_EXECUTABLE(main src/main.cpp ${SOURCES_FILES} ${HEADER_FILES})

add_executable(test src/tests/tests.cpp ${SOURCES_FILES} ${HEADER_FILES})
target_link_libraries(test gtest gtest_main  ${Boost_LIBRARIES})


# set(CMAKE_SHARED_LINKER_FLAGS "-fopenmp -fPIC")

# add_library (bvpodecore OBJECT ${SOURCES_FILES})
# add_library (bvpodestatic STATIC $<TARGET_OBJECTS:bvpodecore>)
# add_library (bvpodeshared SHARED $<TARGET_OBJECTS:bvpodecore>)





