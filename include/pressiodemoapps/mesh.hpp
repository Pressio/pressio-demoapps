/*
//@HEADER
// ************************************************************************
//
// mesh.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef PRESSIODEMOAPPS_MESH_HPP_
#define PRESSIODEMOAPPS_MESH_HPP_

#include <cmath>
#include <limits>
#include <fstream>
#include <array>
#include <iostream>
#include <sstream>

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./impl/mesh_ccu.hpp"

namespace pressiodemoapps{

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class T>
T loadCellCenterUniformMesh(const std::string & meshFilesPath)
{
  return T(meshFilesPath);
}

#else

template<class scalar_type = double>
using cellcentered_uniform_mesh_eigen_type = impl::CellCenteredUniformMesh<
  scalar_type,
  int32_t, // ordinal type
  Eigen::Matrix<scalar_type, -1, 1>, // coordinate storage type
  Eigen::Matrix<int32_t, -1, -1, Eigen::RowMajor> // graph type
  >;

template<class scalar_type = double>
auto load_cellcentered_uniform_mesh_eigen(const std::string & meshFilesPath)
{
  return cellcentered_uniform_mesh_eigen_type<scalar_type>(meshFilesPath);
}

#endif

}//end namespace pressiodemoapps
#endif
