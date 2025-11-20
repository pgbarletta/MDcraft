cmake_minimum_required(VERSION 3.9)
project (mdcraft C CXX)

set(CMAKE_CXX_FLAGS_DEBUG      "-g")
set(CMAKE_CXX_FLAGS_RELEASE   "-O3")
# Select build type (Debug/Release) interactively
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Default build type at CMake configure time")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release")

set(CMAKE_CXX_STANDARD 17 CACHE STRING "Default CXX standard at CMake configure time")
set_property(CACHE CMAKE_CXX_STANDARD PROPERTY STRINGS "11;14;17")

set(mdcraft_VERSION_MAJOR         1)
set(mdcraft_VERSION_MINOR         0)
set(mdcraft_VERSION_PATCH         0)
set(mdcraft_VERSION "${mdcraft_VERSION_MAJOR},${mdcraft_VERSION_MINOR},${mdcraft_VERSION_PATCH}")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;")

option(mdcraft_ENABLE_MPI            "Distributed computing with MPI"                OFF)
option(mdcraft_ENABLE_PYTHON         "Python support"                                OFF)
option(mdcraft_ENABLE_TESTS          "Build tests"                                   OFF)
option(mdcraft_ENABLE_PROBLEMS       "Build problems"                                OFF)
option(mdcraft_ENABLE_BPNN           "Behler-Parrinello NN potentials support"       OFF)
option(mdcraft_ENABLE_DeePMD         "DeePmd NN potentials support"                  OFF)
option(mdcraft_ENABLE_MLIP4          "MTP potentials support from MLIP4"             OFF)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions(-DDEBUG)
endif()

# chosse Multi-threading option
set(mdcraft_THREADS_TYPE "STD" CACHE STRING "Choose thread style: STD or TBB")
set_property(CACHE mdcraft_THREADS_TYPE PROPERTY STRINGS STD TBB)
if (${mdcraft_THREADS_TYPE} STREQUAL TBB)
  message(STATUS "Use TBB")
  set(mdcraft_ENABLE_TBB ON)
else ()
  message(STATUS "Use standard threads")
  set(mdcraft_ENABLE_TBB OFF)
endif ()

find_package(Eigen3)

if(mdcraft_ENABLE_DeePMD)
  find_package(DeePMD REQUIRED)
endif()

if(mdcraft_ENABLE_MLIP4)
  find_package(mlip4 REQUIRED)
endif()

if(mdcraft_ENABLE_PYTHON)
  find_package(pybind11 REQUIRED)
  add_subdirectory(python)
endif()

add_subdirectory(mdcraft)

if(mdcraft_ENABLE_TESTS)
add_subdirectory(test)
endif()

if(mdcraft_ENABLE_PROBLEMS)
add_subdirectory(problems)
endif()
