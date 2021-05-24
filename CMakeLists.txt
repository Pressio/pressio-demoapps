
cmake_minimum_required(VERSION 3.18.0)
project(pressiodemoapps CXX)
set (CMAKE_CXX_STANDARD 11)

#=====================================================================
# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed.
    Please make a new directory (called a build directory) and run CMake from there.
    You may need to remove CMakeCache.txt.")
endif()

#=====================================================================
# check build type
if (NOT CMAKE_BUILD_TYPE)
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
# c++ standard
#=====================================================================
# commands to test if compiler supports standard
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORT_CPP11)
if(NOT COMPILER_SUPPORT_CPP11)
  message(FATAL_ERROR "Compiler does not support -std=c++11. This is required.")
endif()
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#=====================================================================
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/tpls/eigen/eigen3
  ${CMAKE_CURRENT_SOURCE_DIR}/tpls/pressio/include)

#=====================================================================
# tests are enabled if PRESSIODEMOAPPS_ENABLE_TESTS = ON
if(PRESSIODEMOAPPS_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()





#include_directories(
#  ${CMAKE_CURRENT_SOURCE_DIR}/tpls/eigen/eigen3
#  ${CMAKE_CURRENT_SOURCE_DIR}/tpls/pressio/include
#  ${CMAKE_CURRENT_SOURCE_DIR})
#add_executable(reconstructFullFom ${CMAKE_CURRENT_SOURCE_DIR}/src/main_reconstruct_fom_state.cc)
#add_executable(reconstructFullFomOd ${CMAKE_CURRENT_SOURCE_DIR}/src/main_reconstruct_fom_state_od.cc)
#add_executable(fom ${CMAKE_CURRENT_SOURCE_DIR}/src/main_fom.cc)
#add_executable(rom ${CMAKE_CURRENT_SOURCE_DIR}/src/main_rom.cc)
#add_executable(hrrom ${CMAKE_CURRENT_SOURCE_DIR}/src/main_hr_rom.cc)
#add_executable(odrom ${CMAKE_CURRENT_SOURCE_DIR}/src/main_od_rom.cc)
