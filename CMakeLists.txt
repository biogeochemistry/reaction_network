CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(reaction_network)

set(CMAKE_CXX_COMPILER $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/bin/mpicxx)
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/libs/cmake-modules")

FIND_PACKAGE(PETSc REQUIRED)

INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
ADD_DEFINITIONS(${PETSC_DEFINITIONS})

set(EIGEN_DIR "${CMAKE_SOURCE_DIR}/libs/eigen/")
set(GNUPLOT_DIR "${CMAKE_SOURCE_DIR}/libs/gnuplot/")
include_directories(
  ${EIGEN_DIR}
  ${GNUPLOT_DIR}
  ${CMAKE_SOURCE_DIR})


set(HEADER_FILES OdeSolver.hpp Runge-Kutta.hpp BoundaryConditions.hpp SecondOrderOde.hpp FiniteDifferenceGrid.hpp BvpOde.hpp Node.hpp)
set(SOURCES_FILES OdeSolver.cpp Runge-Kutta.cpp BoundaryConditions.cpp BvpOde.cpp FiniteDifferenceGrid.cpp)


ADD_EXECUTABLE(main main.cpp ${SOURCES_FILES} ${HEADER_FILES})
TARGET_LINK_LIBRARIES(main ${PETSC_LIBRARIES})



# CMAKE_MINIMUM_REQUIRED(VERSION 2.8)


# ########################################################################
# #
# # Project-wide settings

# # Name of the project.
# #
# # CMake files in this project can refer to the root source directory
# # as ${gtest_SOURCE_DIR} and to the root binary directory as
# # ${gtest_BINARY_DIR}.
# # Language "C" is required for find_package(Threads).
# project(reaction_network CXX)


# # set(GTEST_DIR "${CMAKE_SOURCE_DIR}/libs/gtest/")
# set(EIGEN_DIR "${CMAKE_SOURCE_DIR}/libs/eigen/")
# set(GNUPLOT_DIR "${CMAKE_SOURCE_DIR}/libs/gnuplot/")
# set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")
# FIND_PACKAGE(PETSc REQUIRED)
# set(CMAKE_CXX_COMPILER ${PETSC_DIR}/${PETSC_ARCH}/bin/mpih)


# # Define helper functions and macros used by Google Test.
# # include(${GTEST_DIR}/cmake/internal_utils.cmake)

# # config_compiler_and_linker()  # Defined in internal_utils.cmake.

# # Where Google Test's .h files can be found.
# include_directories(
#   # ${GTEST_DIR}/include
#   # ${GTEST_DIR}
#   ${EIGEN_DIR}
#   ${GNUPLOT_DIR}
#   # ${PETSC_DIR}
#   # ${PETSC_INCLUDES}
#   ${CMAKE_SOURCE_DIR})
# # Where Google Test's libraries can be found.
# # link_directories(${GTEST_DIR}/src)

# INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
# ADD_DEFINITIONS(${PETSC_DEFINITIONS})

# set(HEADER_FILES OdeSolver.hpp Runge-Kutta.hpp BoundaryConditions.hpp SecondOrderOde.hpp FiniteDifferenceGrid.hpp BvpOde.hpp Node.hpp)
# set(SOURCES_FILES OdeSolver.cpp Runge-Kutta.cpp BoundaryConditions.cpp BvpOde.cpp FiniteDifferenceGrid.cpp)

# # add_definitions("-stdlib=libc++ -std=c++11 -lc++ -g -Wall -Wextra -pthread -Wno-missing-field-initializers")

# ########################################################################
# #
# # Defines the gtest & gtest_main libraries.  User tests should link
# # with one of them.

# # Google Test libraries.  We build them using more strict warnings than what
# # are used for other targets, to ensure that gtest can be compiled by a user
# # aggressive about warnings.
# # cxx_library(gtest "${cxx_strict}" ${GTEST_DIR}/src/gtest-all.cc)
# # cxx_library(gtest_main "${cxx_strict}" ${GTEST_DIR}/src/gtest_main.cc)
# # target_link_libraries(gtest_main gtest ${PETSC_LIBRARIES})

# # cxx_executable(name dir lib srcs...)
# #
# # creates a named target that depends on the given lib and is built
# # from the given source files.  dir/name.cc is implicitly included in
# # the source file list.

# ADD_EXECUTABLE(main main.cpp ${SOURCES_FILES} ${HEADER_FILES})

# add_executable(test tests/test.cpp ${SOURCES_FILES} ${HEADER_FILES})
# target_link_libraries(test gtest gtest_main ${PETSC_LIBRARIES})




