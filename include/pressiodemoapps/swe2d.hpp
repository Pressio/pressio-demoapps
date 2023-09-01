/*
//@HEADER
// ************************************************************************
//
// swe2d.hpp
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

#ifndef PRESSIODEMOAPPS_SWE2D_INC_HPP_
#define PRESSIODEMOAPPS_SWE2D_INC_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./impl/custom_bc_holder.hpp"
#include "./adapter_cpp.hpp"
#include "./adapter_py.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// enums identifying the problems
// ----------------------------------------------------------
enum class Swe2d{
  SlipWall,
  CustomBCs
};

}//end namespace pressiodemoapps

// this include is here because needs visiblity of the enums above
#include "./impl/swe_2d_prob_class.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// create default problem
// ----------------------------------------------------------

#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<implswe2d::EigenApp<mesh_t>>
  >
RetType create_swe2d_problem_default_for_py(const mesh_t & meshObj,
					    Swe2d problemEnum,
					    InviscidFluxReconstruction inviscRecEn)
{
  if (problemEnum == Swe2d::SlipWall){
    return RetType(implswe2d::TagProblemSlipWall{}, meshObj, inviscRecEn,
		   ::pressiodemoapps::InviscidFluxScheme::Rusanov,
		   implswe2d::create_vec_with_default_params<typename RetType::scalar_type>());
  }
  else{
    throw std::runtime_error("2D swe: invalid problem enum");
  }
}

#else

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<implswe2d::EigenApp<mesh_t>>
  >
RetType
create_problem_eigen(const mesh_t & meshObj,
		     Swe2d problemEnum,
		     InviscidFluxReconstruction inviscRecEn,
		     const std::unordered_map<std::string, typename mesh_t::scalar_t> & paramsMap
		       = implswe2d::create_map_with_default_params<typename mesh_t::scalar_t>())
{

  if (problemEnum == Swe2d::SlipWall){
    auto paramsVec = implswe2d::param_unord_map_to_vector(paramsMap);
    return RetType(implswe2d::TagProblemSlipWall{}, meshObj, inviscRecEn,
		   ::pressiodemoapps::InviscidFluxScheme::Rusanov, std::move(paramsVec));
  }
  else{
    throw std::runtime_error("2D swe: invalid problem enum");
  }
}
#endif

#if !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<
  class mesh_t,
  class CustomBCsFunctorLeft,
  class CustomBCsFunctorFront,
  class CustomBCsFunctorRight,
  class CustomBCsFunctorBack,
  class BCFunctorsHolderType = impl::CustomBCsHolder<
    CustomBCsFunctorLeft, CustomBCsFunctorFront,
    CustomBCsFunctorRight, CustomBCsFunctorBack>,
  class RetType = PublicProblemEigenMixinCpp<implswe2d::EigenApp<mesh_t, BCFunctorsHolderType>>
  >
RetType
create_problem_eigen(const mesh_t & meshObj,
		     Swe2d problemEnum,
		     InviscidFluxReconstruction recEnum,
		     CustomBCsFunctorLeft && customBCsLeft,
		     CustomBCsFunctorFront && customBCsFront,
		     CustomBCsFunctorRight && customBCsRight,
		     CustomBCsFunctorBack && customBCsBack,
		     const std::unordered_map<std::string, typename mesh_t::scalar_t> & paramsMap)
{

  if (problemEnum != Swe2d::CustomBCs){
    throw std::runtime_error("Swe2d: custom BCs only supported for Swe2d::CustomBCs");
  }
  if (recEnum != InviscidFluxReconstruction::FirstOrder){
    throw std::runtime_error("Swe2d::CustomBCs: only supports InviscidFluxReconstruction::FirstOrder");
  }

  BCFunctorsHolderType bcFuncs(std::forward<CustomBCsFunctorLeft>(customBCsLeft),
			       std::forward<CustomBCsFunctorFront>(customBCsFront),
			       std::forward<CustomBCsFunctorRight>(customBCsRight),
			       std::forward<CustomBCsFunctorBack>(customBCsBack));
  return RetType(meshObj, problemEnum, recEnum, InviscidFluxScheme::Rusanov,
		 std::move(bcFuncs), implswe2d::param_unord_map_to_vector(paramsMap));
}

template<
  class mesh_t,
  class CustomBCsFunctorLeft,
  class CustomBCsFunctorFront,
  class CustomBCsFunctorRight,
  class CustomBCsFunctorBack>
auto create_problem_eigen(const mesh_t & meshObj,
			  Swe2d problemEnum,
			  InviscidFluxReconstruction recEnum,
			  CustomBCsFunctorLeft && customBCsLeft,
			  CustomBCsFunctorFront && customBCsFront,
			  CustomBCsFunctorRight && customBCsRight,
			  CustomBCsFunctorBack && customBCsBack,
			  [[maybe_unused]] int icFlag = 0 )
{

  auto paramsVec = implswe2d::create_vec_with_default_params<typename mesh_t::scalar_t>();
  auto paramsMap = implswe2d::param_vector_to_unord_map(paramsVec);
  return create_problem_eigen(meshObj, problemEnum, recEnum,
			      std::forward<CustomBCsFunctorLeft>(customBCsLeft),
			      std::forward<CustomBCsFunctorFront>(customBCsFront),
			      std::forward<CustomBCsFunctorRight>(customBCsRight),
			      std::forward<CustomBCsFunctorBack>(customBCsBack),
			      paramsMap);
}
#endif


// ----------------------------------------------------------
// custom coeffs
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<implswe2d::EigenApp<mesh_t>>
  >
RetType
// bindings need unique naming or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  create_slip_wall_swe2d_problem_ov1_for_py
#else
  create_slip_wall_swe_2d_problem_eigen
#endif
(const mesh_t & meshObj,
 InviscidFluxReconstruction inviscRecEn,
 typename mesh_t::scalar_t gravity,
 typename mesh_t::scalar_t coriolis,
 typename mesh_t::scalar_t pulseMagnitude)
{

  auto paramsVec = implswe2d::create_vec_with_default_params<typename RetType::scalar_type>();
  paramsVec[implswe2d::param_string_to_index<>("gravity")] = gravity;
  paramsVec[implswe2d::param_string_to_index<>("coriolis")] = coriolis;
  paramsVec[implswe2d::param_string_to_index<>("pulseMagnitude")] = pulseMagnitude;

  return RetType(implswe2d::TagProblemSlipWall{}, meshObj, inviscRecEn,
		 ::pressiodemoapps::InviscidFluxScheme::Rusanov, std::move(paramsVec));
}

}//end namespace pressiodemoapps
#endif
