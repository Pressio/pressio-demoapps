
cmake_minimum_required(VERSION 3.18.0)
project(pressiodemoapps CXX)

# setting where other cmake commands are
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#=====================================================================
# versioning
#=====================================================================
file(READ "version.txt" PDA_VERSION)
message("pressiodemoapps version = ${PDA_VERSION}")

#=====================================================================
# we need c++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# test if compiler supports standard
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++17" COMPILER_SUPPORT_CPP17)
if(NOT COMPILER_SUPPORT_CPP17)
  message(FATAL_ERROR "Compiler does not support -std=c++17. This is required.")
endif()

# force C++17 for pressio too
set(PRESSIO_ENABLE_CXX17 ON)

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

if(PRESSIODEMOAPPS_ENABLE_OPENMP)
  find_package(OpenMP)
  if(NOT OpenMP_CXX_FOUND)
    message(FATAL_ERROR "You enabled OpenMP but I cannot find it")
  endif()
endif()

if(PRESSIODEMOAPPS_ENABLE_BINDINGS)

  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/tpls/eigen/eigen3)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tpls/pybind11)

  if(PRESSIODEMOAPPS_ENABLE_OPENMP)
    link_libraries(OpenMP::OpenMP_CXX)
  endif()

  set(modulename _pressiodemoappsimpl)
  pybind11_add_module(${modulename} ${CMAKE_CURRENT_SOURCE_DIR}/src_py/main_binder.cc)
  target_compile_definitions(${modulename} PRIVATE MODNAME=${modulename})
  target_compile_definitions(${modulename} PRIVATE PRESSIODEMOAPPS_ENABLE_BINDINGS)
  if(PRESSIODEMOAPPS_ENABLE_OPENMP)
    target_compile_definitions(${modulename} PRIVATE PRESSIODEMOAPPS_ENABLE_OPENMP)
  endif()

else()

  if(PRESSIODEMOAPPS_ENABLE_OPENMP)
    add_compile_definitions(PRESSIODEMOAPPS_ENABLE_OPENMP)
    link_libraries(OpenMP::OpenMP_CXX)
  endif()
  #add_compile_definitions(EIGEN_DONT_PARALLELIZE)

  if(PRESSIODEMOAPPS_ENABLE_TESTS)
    enable_testing()
    add_subdirectory(tests_mesh)
    add_subdirectory(tests_cpp)
  endif()

  if(PRESSIODEMOAPPS_ENABLE_PERF)
    add_subdirectory(tests_perf)
  endif()

endif()
