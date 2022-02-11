
cmake_minimum_required(VERSION 3.16.0)
project(pdaCodeWorkFlow CXX)

#=====================================================================
# we need c++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# test if compiler supports standard
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++14" COMPILER_SUPPORT_CPP14)
if(NOT COMPILER_SUPPORT_CPP14)
  message(FATAL_ERROR "Compiler does not support -std=c++14. This is required.")
endif()

#=====================================================================
# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed.
    Please make a new directory (e.g. build directory) and run CMake from there.")
endif()

# default to release if build type is empty
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "${CMAKE_BUILD_TYPE} is not specified, defaulting to Release.")
  set(CMAKE_BUILD_TYPE "Release")
endif()

# convert cmake build type to lower string
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if( NOT cmake_build_type_tolower STREQUAL "debug"
    AND NOT cmake_build_type_tolower STREQUAL "release")
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\".
    Allowed values are Debug, Release (case-insensitive).")
endif()

#=====================================================================

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/yaml-cpp-yaml-cpp-0.7.0)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/../../include
  ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/eigen/eigen3
  ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/pressio/include
  ${CMAKE_CURRENT_SOURCE_DIR}/yaml-cpp-yaml-cpp-0.7.0/include
  )

find_package(OpenMP)
if(NOT OpenMP_CXX_FOUND)
  message(FATAL_ERROR "You enabled OpenMP but I cannot find it")
endif()


set(exename pdaExeWf)
add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
target_compile_definitions(${exename} PRIVATE PRESSIODEMOAPPS_ENABLE_OPENMP)
target_link_libraries(${exename} PRIVATE yaml-cpp OpenMP::OpenMP_CXX)
