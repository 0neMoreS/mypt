cmake_minimum_required(VERSION 3.20)
project(photon_mapping LANGUAGES CXX)

# set(CMAKE_BUILD_TYPE Debug)

option(BUILD_TESTS "build tests" OFF)

# OpenMP
find_package(OpenMP)
if(NOT OpenMP_CXX_FOUND)
    message(WARNING "failed to find OpenMP")
endif()

# Embree3
find_package(embree 3.12 REQUIRED
    PATHS /home/alex-linux/computer-graphics/mypt/src/external/embree-3.12.2.x86_64.linux/lib/cmake/embree-3.12.2
    NO_DEFAULT_PATH)
if (${embree_FOUND})
  message(STATUS "Found Embree")
else()
  message(FATAL_ERROR "Could not find Embree")
endif()

# externals
add_subdirectory("src/external")

# pm(photon_mapping)
add_library(pm INTERFACE)
target_include_directories(pm INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/src/include")
target_compile_features(pm INTERFACE cxx_std_20)
set_target_properties(pm PROPERTIES CXX_EXTENSIONS OFF)

target_link_libraries(pm INTERFACE OpenMP::OpenMP_CXX)
target_link_libraries(pm INTERFACE tinyobjloader)
target_link_libraries(pm INTERFACE embree)

# examples
add_subdirectory("src")
add_subdirectory("tests")