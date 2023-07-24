/*
//@HEADER
// ************************************************************************
//
// euler2d.hpp
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

#ifndef PRESSIODEMOAPPS_EULER2D_INC_HPP_
#define PRESSIODEMOAPPS_EULER2D_INC_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./euler_compute_energy.hpp"
#include "./adapter_cpp.hpp"
#include "./adapter_py.hpp"
#include "./ghost_relative_locations.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// enums identifying the problems
// ----------------------------------------------------------
enum class Euler2d{
  PeriodicSmooth,
  KelvinHelmholtz,
  SedovFull,
  SedovSymmetry,
  Riemann,
  RiemannCustomBCs,
  NormalShock,
  DoubleMachReflection,
  CrossShock,
  testingonlyneumann
};
}//end namespace pressiodemoapps

// this include is here because needs visiblity of the enums above
#include "./impl/euler_2d_prob_class.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// create a default problem
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impleuler2d::EigenApp<mesh_t>>
  >
RetType
// bindings need unique nameing or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  create_euler_2d_problem_default_for_py
#else
  create_problem_eigen
#endif
(const mesh_t & meshObj,
 Euler2d problemEnum,
 InviscidFluxReconstruction recEnum)
{
  return RetType(meshObj, problemEnum, recEnum,
		 InviscidFluxScheme::Rusanov, 1);
}

#if !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<
  class mesh_t,
  class CustomBCsFunctor,
  class RetType = PublicProblemEigenMixinCpp<impleuler2d::EigenApp<mesh_t, CustomBCsFunctor>>
  >
RetType
create_problem_eigen(const mesh_t & meshObj,
		     Euler2d problemEnum,
		     InviscidFluxReconstruction recEnum,
		     const CustomBCsFunctor & customBCs,
		     int icId = 1)
{

  if (problemEnum != Euler2d::RiemannCustomBCs){
    throw std::runtime_error("custom BCs only supported for Euler2d::RiemannCustomBCs");
  }
  if (recEnum != InviscidFluxReconstruction::FirstOrder){
    throw std::runtime_error("RiemannCustomBCs only supports InviscidFluxReconstruction::FirstOrder");
  }

  return RetType(meshObj, problemEnum, recEnum,
		 InviscidFluxScheme::Rusanov, customBCs, icId);
}
#endif

// ----------------------------------------------------------
// create a default problem with specific initial condition
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impleuler2d::EigenApp<mesh_t>>
  >
RetType
// bindings need unique nameing or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  create_euler_2d_problem_ov1_for_py
#else
  create_problem_eigen
#endif
(const mesh_t & meshObj,
 Euler2d problemEnum,
 InviscidFluxReconstruction recEnum,
 int icId)
{
  return RetType(meshObj, problemEnum, recEnum,
		 InviscidFluxScheme::Rusanov, icId);
}


// this crossshock one should really just be experimental
// and it is not documented on the website
template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impleuler2d::EigenApp<mesh_t>>
  >
RetType
// bindings need unique nameing or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  create_euler_2d_cross_shock_problem_for_py
#else
  create_cross_shock_problem_eigen
#endif
(const mesh_t & meshObj,
 InviscidFluxReconstruction recEnum,
 typename mesh_t::scalar_t density,
 typename mesh_t::scalar_t inletXVel,
 typename mesh_t::scalar_t bottomYVel)
{
  return RetType(impleuler2d::TagCrossShock{},
		 meshObj, recEnum,
		 InviscidFluxScheme::Rusanov, 1,
		 inletXVel, bottomYVel, density);
}

}//end namespace pressiodemoapps
#endif
