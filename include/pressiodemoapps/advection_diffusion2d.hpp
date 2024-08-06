/*
//@HEADER
// ************************************************************************
//
// advection_diffusion2d.hpp
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

#ifndef PRESSIODEMOAPPS_ADVECTION_DIFFUSION_2D_HPP_
#define PRESSIODEMOAPPS_ADVECTION_DIFFUSION_2D_HPP_

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
enum class AdvectionDiffusion2d{
  BurgersPeriodic,
  BurgersOutflow
};

}//end namespace pressiodemoapps

// this include is here because needs visiblity of the enums above
#include "./impl/advection_diffusion_2d_prob_class.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// create default problem
// ----------------------------------------------------------
template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impladvdiff2d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// bindings need unique naming or we get error associated with overloads
create_advecdiffusion2d_problem_default_for_py
#else
create_problem_eigen
#endif
(const mesh_t & meshObj,
 AdvectionDiffusion2d problemEnum,
 InviscidFluxReconstruction inviscidFluxRecEnum,
 ViscousFluxReconstruction viscFluxRecEnum)
{

  using scalar_t = typename mesh_t::scalar_t;
  auto defaultPhsParams = impladvdiff2d::defaultPhysicalParams<scalar_t>;
  const auto defaultICParams = (problemEnum == AdvectionDiffusion2d::BurgersPeriodic || problemEnum == AdvectionDiffusion2d::BurgersOutflow)
    ? impladvdiff2d::defaultInitCondParams<scalar_t> : std::vector<scalar_t>();

  return RetType(meshObj,
  		 problemEnum,
		 inviscidFluxRecEnum,
		 InviscidFluxScheme::Rusanov,
		 viscFluxRecEnum,
		 defaultICParams,
		 defaultPhsParams);
}

// ----------------------------------------------------------
// create problem with custom coefficients
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impladvdiff2d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// bindings need unique naming or we get error associated with overloads
create_burgers_2d_problem_ov1_for_py
#else
create_problem_eigen
#endif
(const mesh_t & meshObj,
 AdvectionDiffusion2d problemEnum,
 InviscidFluxReconstruction inviscidFluxRecEnum,
 ViscousFluxReconstruction viscFluxRecEnum,
 const std::unordered_map<std::string, typename mesh_t::scalar_t> & userParams)
{

  // only one IC for now
  int icFlag = 1;

  using scalar_t = typename mesh_t::scalar_t;
  auto physParamsVec = impladvdiff2d::defaultPhysicalParams<scalar_t>;
  auto icParamsVec   = impladvdiff2d::defaultInitCondParams<scalar_t>;
  impladvdiff2d::replace_params_from_map_if_present(physParamsVec, icParamsVec, problemEnum, icFlag, userParams);

  return RetType(meshObj,
		 problemEnum,
		 inviscidFluxRecEnum,
		 InviscidFluxScheme::Rusanov,
		 viscFluxRecEnum,
		 icParamsVec,
		 physParamsVec);
}

//
// custom BCs
//

#if !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class BCsFuncL, class BCsFuncF, class BCsFuncR, class BCsFuncB>
auto create_problem_eigen(const mesh_t & meshObj,
			  AdvectionDiffusion2d problemEnum,
			  InviscidFluxReconstruction inviscidFluxRecEnum,
			  BCsFuncL && BCsLeft,
			  BCsFuncF && BCsFront,
			  BCsFuncR && BCsRight,
			  BCsFuncB && BCsBack,
			  const int icFlag = 1)
{

  // default viscous flux order
  // TODO: generalize
  ViscousFluxReconstruction viscFluxRecEnum = ViscousFluxReconstruction::FirstOrder;

  using BCFunctorsHolderType = impl::CustomBCsHolder<BCsFuncL, BCsFuncF, BCsFuncR, BCsFuncB>;
  BCFunctorsHolderType bcFuncs(std::forward<BCsFuncL>(BCsLeft),
			       std::forward<BCsFuncF>(BCsFront),
			       std::forward<BCsFuncR>(BCsRight),
			       std::forward<BCsFuncB>(BCsBack));

  using scalar_t = typename mesh_t::scalar_t;
  const auto defaultPhsParams = impladvdiff2d::defaultPhysicalParams<scalar_t>;
  const auto defaultICParams = impladvdiff2d::defaultInitCondParams<scalar_t>;

  using return_type = PublicProblemEigenMixinCpp<impladvdiff2d::EigenApp<mesh_t, BCFunctorsHolderType>>;
  return return_type(meshObj,
		     problemEnum,
		     inviscidFluxRecEnum,
		     InviscidFluxScheme::Rusanov,
		     viscFluxRecEnum,
		     std::move(bcFuncs),
		     defaultPhsParams,
		     defaultICParams);

}
#endif

#if !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class BCsFuncL, class BCsFuncF, class BCsFuncR, class BCsFuncB>
auto create_problem_eigen(const mesh_t & meshObj,
			  AdvectionDiffusion2d problemEnum,
			  InviscidFluxReconstruction inviscidFluxRecEnum,
			  BCsFuncL && BCsLeft,
			  BCsFuncF && BCsFront,
			  BCsFuncR && BCsRight,
			  BCsFuncB && BCsBack,
			  const int icFlag,
			  const std::unordered_map<std::string, typename mesh_t::scalar_t> & userParams)
{

  // default viscous flux order
  // TODO: generalize
  ViscousFluxReconstruction viscFluxRecEnum = ViscousFluxReconstruction::FirstOrder;

  using BCFunctorsHolderType = impl::CustomBCsHolder<BCsFuncL, BCsFuncF, BCsFuncR, BCsFuncB>;
  BCFunctorsHolderType bcFuncs(std::forward<BCsFuncL>(BCsLeft),
			       std::forward<BCsFuncF>(BCsFront),
			       std::forward<BCsFuncR>(BCsRight),
			       std::forward<BCsFuncB>(BCsBack));

  using scalar_t = typename mesh_t::scalar_t;
  auto physParamsVec = impladvdiff2d::defaultPhysicalParams<scalar_t>;
  auto icParamsVec   = impladvdiff2d::defaultInitCondParams<scalar_t>;
  impladvdiff2d::replace_params_from_map_if_present(physParamsVec, icParamsVec, problemEnum, icFlag, userParams);

  using return_type = PublicProblemEigenMixinCpp<impladvdiff2d::EigenApp<mesh_t, BCFunctorsHolderType>>;
  return return_type(meshObj,
		     problemEnum,
		     inviscidFluxRecEnum,
		     InviscidFluxScheme::Rusanov,
		     viscFluxRecEnum,
		     std::move(bcFuncs),
		     icParamsVec,
		     physParamsVec);
}
#endif

} //end namespace pressiodemoapps
#endif
