# =========================================================
# CFD Basics

cmake_minimum_required(VERSION 3.25)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "cmake")

set(ProjectName "cfd-basics")
set(ProjectNamespace "cfd_basics")

project("${ProjectName}")

set(CMAKE_CXX_STANDARD 20)

add_compile_options(
        -Wall -Werror -Wno-error=uninitialized
        $<$<CONFIG:RELEASE>:-O3>
        $<$<CONFIG:RELEASE>:-DNDEBUG>
        $<$<CONFIG:RELWITHDEBINFO>:-O2>
        $<$<CONFIG:RELWITHDEBINFO>:-DNDEBUG>
        $<$<CONFIG:RELWITHDEBINFO>:-g>
        $<$<CONFIG:MINSIZEREL>:-Oz>
        $<$<CONFIG:MINSIZEREL>:-DNDEBUG>
        $<$<CONFIG:DEBUG>:-Og>
        $<$<CONFIG:DEBUG>:-g>
)

# Add module directory and modules
include(FetchContent)

# =========================================================
# gtest

FetchContent_Declare(
        googletest
        URL
        https://github.com/google/googletest/archive/refs/tags/v1.13.0.zip
)

set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(GTEST_FORCE_SHARED_CRT ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

set_property(TARGET gmock PROPERTY FOLDER "gtest")
set_property(TARGET gmock_main PROPERTY FOLDER "gtest")
set_property(TARGET gtest PROPERTY FOLDER "gtest")
set_property(TARGET gtest_main PROPERTY FOLDER "gtest")

# =========================================================
# Eigen

FetchContent_Declare(
        Eigen
        URL
        https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
)

set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(EIGEN_FORCE_SHARED_CRT ON CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DOC OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)
FetchContent_MakeAvailable(Eigen)

# =========================================================
# Third party libraries
include(third_party/nuenv.cmake)

# =========================================================
# Targets

# --- project

add_library("${ProjectName}" STATIC
        src/builder/boundary_field.hpp
        src/builder/cfd_utils.hpp
        src/builder/diffusion.hpp
        src/builder/system.hpp
        src/linalg/gauss_seidel.hpp
        src/linalg/tridiag.hpp
        src/mesh/cell.hpp
        src/mesh/face.hpp
        src/mesh/poly_boundary.hpp
        src/mesh/poly_mesh.hpp
        src/numcal/integration.hpp
        #        src/Core.h
        src/scripts.hpp
        src/solution.hpp
)

target_include_directories("${ProjectName}" PUBLIC
        "${CMAKE_SOURCE_DIR}/src"
)

target_link_libraries("${ProjectName}" PUBLIC
        Eigen3::Eigen
        Nuenv::nuenv
)

target_compile_definitions("${ProjectName}" PRIVATE PROJECT_NAMESPACE=${ProjectNamespace})
set_target_properties("${ProjectName}" PROPERTIES LINKER_LANGUAGE CXX)

# --- main

add_executable("${ProjectName}-main" src/main.cpp)
target_link_libraries("${ProjectName}-main" PUBLIC "${ProjectName}")

# --- tc1

add_executable("${ProjectName}-tc1" tc1/tc1.cpp)
target_link_libraries("${ProjectName}-tc1" PUBLIC "${ProjectName}")

# --- gtest

enable_testing()

add_executable("${ProjectName}-test"
        test/sanity_test.cpp
)

target_link_libraries("${ProjectName}-test" gtest_main "${ProjectName}")

include(GoogleTest)
gtest_discover_tests("${ProjectName}-test")

# =========================================================
