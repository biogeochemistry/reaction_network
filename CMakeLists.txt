########################################################################
# CMake build script for Google Test.
#
# To run the tests for Google Test itself on Linux, use 'make test' or
# ctest.  You can select which tests to run using 'ctest -R regex'.
# For more options, run 'ctest --help'.

# BUILD_SHARED_LIBS is a standard CMake variable, but we declare it here to
# make it prominent in the GUI.
option(BUILD_SHARED_LIBS "Build shared libraries (DLLs)." OFF)

# When other libraries are using a shared version of runtime libraries,
# Google Test also has to use one.
option(  gtest_force_shared_crt
  "Use shared (DLL) run-time lib even when Google Test is built as static lib."
  OFF)

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


# Change this
set(reaction_network ${CMAKE_SOURCE_DIR})


set(GTEST_DIR "${reaction_network}/libs/gtest/")
set(EIGEN_DIR "${reaction_network}/libs/eigen/")
set(GNUPLOT_DIR "${reaction_network}/libs/gnuplot/")


# Define helper functions and macros used by Google Test.
include(${GTEST_DIR}/cmake/internal_utils.cmake)

config_compiler_and_linker()  # Defined in internal_utils.cmake.

# Where Google Test's .h files can be found.
include_directories(
  ${GTEST_DIR}/include
  ${GTEST_DIR}
  ${EIGEN_DIR}
  ${GNUPLOT_DIR}
  ${reaction_network})
# Where Google Test's libraries can be found.
link_directories(${GTEST_DIR}/src)

set(HEADER_FILES OdeSolver.hpp Runge-Kutta.hpp BoundaryConditions.hpp SecondOrderOde.hpp FiniteDifferenceGrid.hpp BvpOde.hpp Node.hpp)
set(SOURCES_FILES OdeSolver.cpp Runge-Kutta.cpp BoundaryConditions.cpp BvpOde.cpp FiniteDifferenceGrid.cpp)

add_definitions("-stdlib=libc++ -std=c++11 -lc++ -g -Wall -Wextra -pthread -Wno-missing-field-initializers")

########################################################################
#
# Defines the gtest & gtest_main libraries.  User tests should link
# with one of them.

# Google Test libraries.  We build them using more strict warnings than what
# are used for other targets, to ensure that gtest can be compiled by a user
# aggressive about warnings.
cxx_library(gtest "${cxx_strict}" ${GTEST_DIR}/src/gtest-all.cc)
cxx_library(gtest_main "${cxx_strict}" ${GTEST_DIR}/src/gtest_main.cc)
target_link_libraries(gtest_main gtest)

# cxx_executable(name dir lib srcs...)
#
# creates a named target that depends on the given lib and is built
# from the given source files.  dir/name.cc is implicitly included in
# the source file list.

ADD_EXECUTABLE(main main.cpp ${SOURCES_FILES} ${HEADER_FILES})

add_executable(test tests/test.cpp ${SOURCES_FILES} ${HEADER_FILES})
target_link_libraries(test gtest gtest_main)




