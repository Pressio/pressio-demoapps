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

/*
create_problem_eigen(meshObj, Swe2d en, InviscidFluxReconstruction);

  - default init condition, default coefficients

create_problem_eigen(meshObj, Swe2d en, InviscidFluxReconstruction, IC = 1);

  - IC: integer specifying the init condition
    - IC = 1: single gaussian pulse, IC = 2: two gaussian pulses
  - default parameters

create_problem_eigen(meshObj, Swe2d en, InviscidFluxReconstruction, IC, params);

  - IC: integer specifying the init condition
    - IC = 1: single gaussian pulse
    - IC = 2: two gaussian pulses
  - params: can contain both physical and params customizing the IC
    - if IC = 1, user can set any of {gravity, coriolis, mag and loc of pulse}
    - if IC = 2, user can set any of {gravity, coriolis, mag and loc for both pulses}
    if a param is not set, we use default value

create_problem_eigen(meshObj, Swe2d en, InviscidFluxReconstruction, customBCs, IC = 0);

  - IC: integer specifying the init condition
    - IC = 1: single gaussian pulse, IC = 2: two gaussian pulses
  - default parameters
  - customBCs

create_problem_eigen(meshObj, Swe2d en, InviscidFluxReconstruction, customBCs, IC, params);

  - IC: integer specifying the init condition
    - IC = 1: single gaussian pulse, IC = 2: two gaussian pulses
  - params: can contain both physical and params customizing the IC
    - if IC = 1, user can set any of {gravity, coriolis, mag and loc of pulse}
    - if IC = 2, user can set any of {gravity, coriolis, mag and loc for both pulses}
    if a param is not set, we use default value
  - customBCs
*/


// ----------------------------------------------------------
// default problem:
//	default init condition, default phys coefficients
// ----------------------------------------------------------

#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class MeshType, class ReturnType>
ReturnType create_swe2d_problem_default_for_py(const MeshType & meshObj,
					    Swe2d problemEnum,
					    InviscidFluxReconstruction inviscRecEn)
{

  if (problemEnum != Swe2d::SlipWall){
    throw std::runtime_error("2D swe: create_swe2d_problem_default_for_py: only supports Swed2d::SlipWall");
  }

  using sc_t = typename MeshType::scalar_t;
  const int icFlag = 1;
  return ReturnType(implswe2d::TagProblemSlipWall{}, meshObj, inviscRecEn,
		    InviscidFluxScheme::Rusanov, icFlag,
		    implswe2d::defaultInitCondParams<sc_t>,
		    implswe2d::defaultPhysicalParams<sc_t>);
}
#else

template<class MeshType>
auto create_problem_eigen(const MeshType & meshObj,
			  Swe2d problemEnum,
			  InviscidFluxReconstruction inviscRecEn,
			  const int icFlag = 1)
{
  // preconditions
  if (problemEnum != Swe2d::SlipWall){
    throw std::runtime_error("2D swe: create_problem_eigen: overload only supporting Swed2d::SlipWall");
  }
  if (!implswe2d::valid_ic_flag<>(icFlag)){
    throw std::runtime_error("2D swe: create_problem_eigen: for Swed2d::SlipWall, invalid icFlag");
  }

  // create problem
  using sc_t = typename MeshType::scalar_t;
  using return_type = PublicProblemEigenMixinCpp<implswe2d::EigenApp<MeshType>>;
  return return_type(implswe2d::TagProblemSlipWall{}, meshObj, inviscRecEn,
		     InviscidFluxScheme::Rusanov, icFlag,
		     implswe2d::defaultInitCondParams<sc_t>,
		     implswe2d::defaultPhysicalParams<sc_t>);
}

template<class MeshType>
auto create_problem_eigen(const MeshType & meshObj,
			  Swe2d problemEnum,
			  InviscidFluxReconstruction inviscRecEn,
			  const int icFlag,
			  const std::unordered_map<std::string, typename MeshType::scalar_t> & userParams)
{
  // preconditions
  if (problemEnum != Swe2d::SlipWall){
    throw std::runtime_error("2D swe: create_problem_eigen: overload only supporting Swed2d::SlipWall");
  }
  if (!implswe2d::valid_ic_flag<>(icFlag)){
    throw std::runtime_error("2D swe: create_problem_eigen: for Swed2d::SlipWall, invalid icFlag");
  }

  if (!implswe2d::contains_valid_parameter_names(userParams)){
    throw std::runtime_error("2D swe: create_problem_eigen: one or more params in user-provided map is invalid");
  }

  using sc_t = typename MeshType::scalar_t;
  using return_type = PublicProblemEigenMixinCpp<implswe2d::EigenApp<MeshType>>;

  // create vec with phys and IC params and then replace those with user defined ones
  auto physParamsVec = implswe2d::defaultPhysicalParams<sc_t>;
  auto icParamsVec   = implswe2d::defaultInitCondParams<sc_t>;
  implswe2d::replace_params_from_map_if_present(physParamsVec, icParamsVec, icFlag, userParams);

  return return_type(implswe2d::TagProblemSlipWall{}, meshObj, inviscRecEn,
		     InviscidFluxScheme::Rusanov, icFlag, icParamsVec, physParamsVec);
}
#endif

#if !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class MeshType, class BCsFuncL, class BCsFuncF, class BCsFuncR, class BCsFuncB>
auto create_problem_eigen(const MeshType & meshObj,
			  Swe2d problemEnum,
			  InviscidFluxReconstruction recEnum,
			  BCsFuncL && BCsLeft,
			  BCsFuncF && BCsFront,
			  BCsFuncR && BCsRight,
			  BCsFuncB && BCsBack,
			  const int icFlag = 1)
{

  if (problemEnum != Swe2d::CustomBCs){
    throw std::runtime_error("2D swe: passing custom BCs requires Swed2d::CustomBCs");
  }
  if (!implswe2d::valid_ic_flag<>(icFlag)){
    throw std::runtime_error("2D swe: for Swed2d::CustomBCs, invalid icFlag");
  }

  using sc_t = typename MeshType::scalar_t;
  using BCFunctorsHolderType = impl::CustomBCsHolder<BCsFuncL, BCsFuncF, BCsFuncR, BCsFuncB>;
  using return_type = PublicProblemEigenMixinCpp<implswe2d::EigenApp<MeshType, BCFunctorsHolderType>>;

  BCFunctorsHolderType bcFuncs(std::forward<BCsFuncL>(BCsLeft),
			       std::forward<BCsFuncF>(BCsFront),
			       std::forward<BCsFuncR>(BCsRight),
			       std::forward<BCsFuncB>(BCsBack));

  return return_type(implswe2d::TagProblemCustomBCs{}, meshObj, recEnum,
		     InviscidFluxScheme::Rusanov, std::move(bcFuncs), icFlag,
		     implswe2d::defaultInitCondParams<sc_t>,
		     implswe2d::defaultPhysicalParams<sc_t>);
}

template<class MeshType, class BCsFuncL, class BCsFuncF, class BCsFuncR, class BCsFuncB>
auto create_problem_eigen(const MeshType & meshObj,
			  Swe2d problemEnum,
			  InviscidFluxReconstruction recEnum,
			  BCsFuncL && BCsLeft,
			  BCsFuncF && BCsFront,
			  BCsFuncR && BCsRight,
			  BCsFuncB && BCsBack,
			  const int icFlag,
			  const std::unordered_map<std::string, typename MeshType::scalar_t> & userParams)
{

  // preconditions
  if (problemEnum != Swe2d::CustomBCs){
    throw std::runtime_error("2D swe: passing custom BCs requires Swed2d::CustomBCs");
  }
  if (!implswe2d::valid_ic_flag<>(icFlag)){
    throw std::runtime_error("2D swe: for Swed2d::CustomBCs, invalid icFlag");
  }
  if (!implswe2d::contains_valid_parameter_names(userParams)){
    throw std::runtime_error("2D swe: create_problem_eigen: one or more params in user-provided map is invalid");
  }

  using sc_t = typename MeshType::scalar_t;
  using BCFunctorsHolderType = impl::CustomBCsHolder<BCsFuncL, BCsFuncF, BCsFuncR, BCsFuncB>;
  using return_type = PublicProblemEigenMixinCpp<implswe2d::EigenApp<MeshType, BCFunctorsHolderType>>;

  auto physParamsVec = implswe2d::defaultPhysicalParams<sc_t>;
  auto icParamsVec   = implswe2d::defaultInitCondParams<sc_t>;
  implswe2d::replace_params_from_map_if_present(physParamsVec, icParamsVec, icFlag, userParams);

  BCFunctorsHolderType bcFuncs(std::forward<BCsFuncL>(BCsLeft),
			       std::forward<BCsFuncF>(BCsFront),
			       std::forward<BCsFuncR>(BCsRight),
			       std::forward<BCsFuncB>(BCsBack));

  return return_type(implswe2d::TagProblemCustomBCs{}, meshObj, recEnum,
		     InviscidFluxScheme::Rusanov, std::move(bcFuncs),
		     icFlag, icParamsVec, physParamsVec);
}
#endif

// ----------------------------------------------------------
// custom coeffs: legacy one to deprecate
// ----------------------------------------------------------

template<
  class MeshType,
  class ReturnType = PublicProblemEigenMixinCpp<implswe2d::EigenApp<MeshType>>
  >
ReturnType
// bindings need unique naming or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  create_slip_wall_swe2d_problem_ov1_for_py
#else
  create_slip_wall_swe_2d_problem_eigen
#endif
(const MeshType & meshObj,
 InviscidFluxReconstruction inviscRecEn,
 typename MeshType::scalar_t gravity,
 typename MeshType::scalar_t coriolis,
 typename MeshType::scalar_t pulseMagnitude)
{

  const int icFlag = 1;

  using sc_t = typename MeshType::scalar_t;
  auto physParamsVec = implswe2d::defaultPhysicalParams<sc_t>;
  auto icParamsVec   = implswe2d::defaultInitCondParams<sc_t>;
  physParamsVec[implswe2d::phys_param_string_to_index<>("gravity")] = gravity;
  physParamsVec[implswe2d::phys_param_string_to_index<>("coriolis")] = coriolis;
  icParamsVec[implswe2d::ic_param_string_to_index<>(icFlag, "pulseMagnitude")] = pulseMagnitude;

  return ReturnType(implswe2d::TagProblemSlipWall{}, meshObj, inviscRecEn,
		    ::pressiodemoapps::InviscidFluxScheme::Rusanov,
		    icFlag, icParamsVec, physParamsVec);
}

}//end namespace pressiodemoapps
#endif
