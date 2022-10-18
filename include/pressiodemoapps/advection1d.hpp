/*
//@HEADER
// ************************************************************************
//
// advection1d.hpp
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

#ifndef PRESSIODEMOAPPS_ADVECTION_1D_HPP_
#define PRESSIODEMOAPPS_ADVECTION_1D_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./adapter_cpp.hpp"
#include "./adapter_py.hpp"
#include "mypi.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// enums identifying the problems
// ----------------------------------------------------------
enum class Advection1d{
  PeriodicLinear
};

}//end namespace pressiodemoapps

// this include is here because needs visiblity of the enums above
#include "./impl/advection_1d_prob_class.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// default problem
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impladvection1d::EigenApp<mesh_t>>
  >
RetType
// bindings need unique naming or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
create_linear_advection_1d_problem_default_for_py
#else
create_problem_eigen
#endif
(const mesh_t & meshObj,
 Advection1d problemEnum,
 InviscidFluxReconstruction inviscidFluxRecEnum)
{

  using scalar_t = typename mesh_t::scalar_t;
  constexpr auto defaultVel = static_cast<scalar_t>(1);
  constexpr int defaultIc = 1;
  if (problemEnum == ::pressiodemoapps::Advection1d::PeriodicLinear){
    return RetType(impladvection1d::TagLinearAdvection{},
		   meshObj, inviscidFluxRecEnum,
		   InviscidFluxScheme::Rusanov,
		   defaultVel, defaultIc);
  }
  else{
    throw std::runtime_error("advection: invalid problem enum");
  }
}

// ----------------------------------------------------------
// problem with custom velocity
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impladvection1d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// bindings need unique naming or we get error associated with overloads
create_linear_advection_1d_problem_ov1_for_py
#else
create_linear_advection_1d_problem_eigen
#endif
(const mesh_t & meshObj,
 InviscidFluxReconstruction inviscidFluxRecEnum,
 InviscidFluxScheme inviscidFluxScheme,
 typename mesh_t::scalar_t velocity)
{

  constexpr int defaultIc = 1;
  return RetType(impladvection1d::TagLinearAdvection{},
		 meshObj, inviscidFluxRecEnum,
		 inviscidFluxScheme, velocity, defaultIc);
}

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impladvection1d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// bindings need unique naming or we get error associated with overloads
create_linear_advection_1d_problem_ov2_for_py
#else
create_linear_advection_1d_problem_eigen
#endif
(const mesh_t & meshObj,
 InviscidFluxReconstruction inviscidFluxRecEnum,
 typename mesh_t::scalar_t velocity,
 int ic=1)
{
  return RetType(impladvection1d::TagLinearAdvection{},
		 meshObj, inviscidFluxRecEnum,
		 InviscidFluxScheme::Rusanov,
		 velocity, ic);
}

}//end namespace pressiodemoapps
#endif
