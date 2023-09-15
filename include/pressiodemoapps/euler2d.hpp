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
#include "./impl/custom_bc_holder.hpp"
#include "./euler_compute_energy.hpp"
#include "./adapter_cpp.hpp"
#include "./adapter_py.hpp"

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
  class MeshType,
  class RetType = PublicProblemEigenMixinCpp<impleuler2d::EigenApp<MeshType>>
  >
RetType
// bindings need unique nameing or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  create_euler_2d_problem_default_for_py
#else
  create_problem_eigen
#endif
(const MeshType & meshObj,
 Euler2d problemEnum,
 InviscidFluxReconstruction recEnum)
{
  // currently only normal shock and Riemann allow user IC parametrization
  // while the other problem are hard-wired
  using sc_t = typename MeshType::scalar_t;
  const auto defaultICParams = (problemEnum == Euler2d::NormalShock || problemEnum == Euler2d::Riemann)
    ? impleuler2d::defaultInitCondParams<sc_t> : std::vector<sc_t>();

  return RetType(meshObj, problemEnum, recEnum,
		 InviscidFluxScheme::Rusanov, 1,
		 defaultICParams,
		 impleuler2d::defaultPhysicalParams<sc_t>);
}

// ----------------------------------------------------------
// create a default problem with specific initial condition
// ----------------------------------------------------------

template<
  class MeshType,
  class RetType = PublicProblemEigenMixinCpp<impleuler2d::EigenApp<MeshType>>
  >
RetType
// bindings need unique nameing or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  create_euler_2d_problem_ov1_for_py
#else
  create_problem_eigen
#endif
(const MeshType & meshObj,
 Euler2d problemEnum,
 InviscidFluxReconstruction recEnum,
 int icFlag)
{

  // currently only normal shock and Riemann allow user IC parametrization
  // while the other problem are hard-wired
  using sc_t = typename MeshType::scalar_t;
  const auto defaultICParams = (problemEnum == Euler2d::NormalShock || problemEnum == Euler2d::Riemann)
    ? impleuler2d::defaultInitCondParams<sc_t> : std::vector<sc_t>();

  return RetType(meshObj, problemEnum, recEnum,
		 InviscidFluxScheme::Rusanov, icFlag,
		 defaultICParams,
		 impleuler2d::defaultPhysicalParams<sc_t>);
}

#if !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class MeshType>
auto create_problem_eigen(const MeshType & meshObj,
			  Euler2d problemEnum,
			  InviscidFluxReconstruction recEnum,
			  const int icFlag,
			  const std::unordered_map<std::string, typename MeshType::scalar_t> & userParams)
{

  if (problemEnum != Euler2d::Riemann && problemEnum != Euler2d::NormalShock){
    throw std::runtime_error("Euler2d: custom parametrization only valid for Euler2d::{Riemann, NormalShock}");
  }
  if (!impleuler2d::valid_ic_flag<>(problemEnum, icFlag)){
    throw std::runtime_error("Euler2d: invalid icFlag for the given problem enum");
  }

  using sc_t = typename MeshType::scalar_t;
  using return_type = PublicProblemEigenMixinCpp<impleuler2d::EigenApp<MeshType>>;

  auto physParamsVec = impleuler2d::defaultPhysicalParams<sc_t>;
  auto icParamsVec   = impleuler2d::defaultInitCondParams<sc_t>;
  impleuler2d::replace_params_from_map_if_present(physParamsVec, icParamsVec, problemEnum, icFlag, userParams);

  return return_type(meshObj, problemEnum, recEnum, InviscidFluxScheme::Rusanov,
		     icFlag, icParamsVec, physParamsVec);
}
#endif


//
// custom BCs
//
#if !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class MeshType, class BCsFuncL, class BCsFuncF, class BCsFuncR, class BCsFuncB>
auto create_problem_eigen(const MeshType & meshObj,
			  Euler2d problemEnum,
			  InviscidFluxReconstruction recEnum,
			  BCsFuncL && BCsLeft,
			  BCsFuncF && BCsFront,
			  BCsFuncR && BCsRight,
			  BCsFuncB && BCsBack,
			  const int icFlag = 1)
{

  if (problemEnum != Euler2d::Riemann && problemEnum != Euler2d::NormalShock){
    throw std::runtime_error("Euler2d: custom BCs only valid for Euler2d::{Riemann, NormalShock}");
  }
  if (!impleuler2d::valid_ic_flag<>(problemEnum, icFlag)){
    throw std::runtime_error("Euler2d: invalid icFlag for the given problem enum");
  }
  if (recEnum != InviscidFluxReconstruction::FirstOrder){
    throw std::runtime_error("Euler2d: using custom BCs is only supported with InviscidFluxReconstruction::FirstOrder");
  }

  using sc_t = typename MeshType::scalar_t;
  using BCFunctorsHolderType = impl::CustomBCsHolder<BCsFuncL, BCsFuncF, BCsFuncR, BCsFuncB>;
  using return_type = PublicProblemEigenMixinCpp<impleuler2d::EigenApp<MeshType, BCFunctorsHolderType>>;

  const auto defaultICParams = impleuler2d::defaultInitCondParams<sc_t>;
  BCFunctorsHolderType bcFuncs(std::forward<BCsFuncL>(BCsLeft),
			       std::forward<BCsFuncF>(BCsFront),
			       std::forward<BCsFuncR>(BCsRight),
			       std::forward<BCsFuncB>(BCsBack));
  return return_type(meshObj, problemEnum, recEnum, InviscidFluxScheme::Rusanov,
		     std::move(bcFuncs), icFlag,
		     defaultICParams,
		     impleuler2d::defaultPhysicalParams<sc_t>);
}

template<class MeshType, class BCsFuncL, class BCsFuncF, class BCsFuncR, class BCsFuncB>
auto create_problem_eigen(const MeshType & meshObj,
			  Euler2d problemEnum,
			  InviscidFluxReconstruction recEnum,
			  BCsFuncL && BCsLeft,
			  BCsFuncF && BCsFront,
			  BCsFuncR && BCsRight,
			  BCsFuncB && BCsBack,
			  const int icFlag,
			  const std::unordered_map<std::string, typename MeshType::scalar_t> & userParams)
{

  if (problemEnum != Euler2d::Riemann && problemEnum != Euler2d::NormalShock){
    throw std::runtime_error("Euler2d: custom BCs only valid for Euler2d::{Riemann, NormalShock}");
  }
  if (!impleuler2d::valid_ic_flag<>(problemEnum, icFlag)){
    throw std::runtime_error("Euler2d: invalid icFlag for the given problem enum");
  }
  if (recEnum != InviscidFluxReconstruction::FirstOrder){
    throw std::runtime_error("Euler2d: using custom BCs is only supported with InviscidFluxReconstruction::FirstOrder");
  }

  using sc_t = typename MeshType::scalar_t;
  using BCFunctorsHolderType = impl::CustomBCsHolder<BCsFuncL, BCsFuncF, BCsFuncR, BCsFuncB>;
  using return_type = PublicProblemEigenMixinCpp<impleuler2d::EigenApp<MeshType, BCFunctorsHolderType>>;

  auto physParamsVec = impleuler2d::defaultPhysicalParams<sc_t>;
  auto icParamsVec   = impleuler2d::defaultInitCondParams<sc_t>;
  impleuler2d::replace_params_from_map_if_present(physParamsVec, icParamsVec, problemEnum, icFlag, userParams);

  BCFunctorsHolderType bcFuncs(std::forward<BCsFuncL>(BCsLeft),
			       std::forward<BCsFuncF>(BCsFront),
			       std::forward<BCsFuncR>(BCsRight),
			       std::forward<BCsFuncB>(BCsBack));
  return return_type(meshObj, problemEnum, recEnum,
		     InviscidFluxScheme::Rusanov, std::move(bcFuncs),
		     icFlag, icParamsVec, physParamsVec);
}
#endif


// this crossshock one should really just be experimental
// and it is not documented on the website
template<
  class MeshType,
  class RetType = PublicProblemEigenMixinCpp<impleuler2d::EigenApp<MeshType>>
  >
RetType
// bindings need unique nameing or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  create_euler_2d_cross_shock_problem_for_py
#else
  create_cross_shock_problem_eigen
#endif
(const MeshType & meshObj,
 InviscidFluxReconstruction recEnum,
 typename MeshType::scalar_t density,
 typename MeshType::scalar_t inletXVel,
 typename MeshType::scalar_t bottomYVel)
{

  using sc_t = typename MeshType::scalar_t;
  auto physParamsVec = impleuler2d::defaultPhysicalParams<sc_t>;
  auto icParamsVec   = impleuler2d::defaultInitCondParams<sc_t>;

  const auto i1 = impleuler2d::ic_param_string_to_index<>(Euler2d::CrossShock, 1, "crossShockDensity");
  const auto i2 = impleuler2d::ic_param_string_to_index<>(Euler2d::CrossShock, 1, "crossShockInletXVel");
  const auto i3 = impleuler2d::ic_param_string_to_index<>(Euler2d::CrossShock, 1, "crossShockBottomYVel");
  icParamsVec[i1] = density;
  icParamsVec[i2] = inletXVel;
  icParamsVec[i3] = bottomYVel;

  return RetType(impleuler2d::TagCrossShock{}, meshObj, recEnum,
		 InviscidFluxScheme::Rusanov, 1, icParamsVec, physParamsVec);
}
#if !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class MeshType>
auto create_cross_shock_problem_eigen(const MeshType & meshObj,
				       InviscidFluxReconstruction recEnum)
{
  using sc_t = typename MeshType::scalar_t;
  using result_type = PublicProblemEigenMixinCpp<impleuler2d::EigenApp<MeshType>>;
  return result_type(impleuler2d::TagCrossShock{}, meshObj, recEnum,
		     InviscidFluxScheme::Rusanov, 1,
		     impleuler2d::defaultInitCondParams<sc_t>,
		     impleuler2d::defaultPhysicalParams<sc_t>);
}
#endif

}//end namespace pressiodemoapps
#endif
