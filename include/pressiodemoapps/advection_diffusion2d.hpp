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
  if (problemEnum == ::pressiodemoapps::AdvectionDiffusion2d::BurgersPeriodic ||
      problemEnum == ::pressiodemoapps::AdvectionDiffusion2d::BurgersOutflow)
  {
    // default parameters
    const auto icPulseMagnitude = static_cast<scalar_t>(0.5);
    const auto icSpread = static_cast<scalar_t>(0.15);
    const auto diffusion = static_cast<scalar_t>(0.00001);
    const auto icCenterX = static_cast<scalar_t>(0.0);
    const auto icCenterY = static_cast<scalar_t>(-0.2);

    return RetType(meshObj,
		   problemEnum,
		   inviscidFluxRecEnum,
		   InviscidFluxScheme::Rusanov,
		   viscFluxRecEnum,
		   icPulseMagnitude, icSpread, diffusion,
		   icCenterX, icCenterY);
  }
  else{
    throw std::runtime_error("advection-diffusion2d: invalid problem enum");
  }
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
create_periodic_burgers_2d_problem_ov1_for_py
#else
create_periodic_burgers_2d_problem_eigen
#endif
(const mesh_t & meshObj,
 AdvectionDiffusion2d problemEnum,
 InviscidFluxReconstruction inviscidFluxRecEnum,
 ViscousFluxReconstruction viscFluxRecEnum,
 typename mesh_t::scalar_t icPulseMagnitude,
 typename mesh_t::scalar_t icSpread,
 typename mesh_t::scalar_t diffusion,
 typename mesh_t::scalar_t icCenterX,
 typename mesh_t::scalar_t icCenterY)
{

  return RetType(meshObj,
		 problemEnum,
		 inviscidFluxRecEnum,
		 InviscidFluxScheme::Rusanov,
		 viscFluxRecEnum,
		 icPulseMagnitude, icSpread, diffusion,
		 icCenterX, icCenterY);
}

//
// custom BCs
//

// #if !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// template<class mesh_t, class BCsFuncL, class BCsFuncF, class BCsFuncR, class BCsFuncB>
// auto create_problem_eigen
// (const mesh_t & meshObj,
//  AdvectionDiffusion2d problemEnum,
//  InviscidFluxReconstruction inviscidFluxRecEnum,
//  ViscousFluxReconstruction viscFluxRecEnum,
//  BCsFuncL && BCsLeft,
//  BCsFuncF && BCsFront,
//  BCsFuncR && BCsRight,
//  BCsFuncB && BCsBack)
// {

//   using scalar_t = typename mesh_t::scalar_t;
//   using BCFunctorsHolderType = impl::CustomBCsHolder<BCsFuncL, BCsFuncF, BCsFuncR, BCsFuncB>;
//   BCFunctorsHolderType bcFuncs(std::forward<BCsFuncL>(BCsLeft),
// 			       std::forward<BCsFuncF>(BCsFront),
// 			       std::forward<BCsFuncR>(BCsRight),
// 			       std::forward<BCsFuncB>(BCsBack));

//   // default parameters
//   const auto icPulseMagnitude = static_cast<scalar_t>(0.5);
//   const auto icSpread = static_cast<scalar_t>(0.15);
//   const auto diffusion = static_cast<scalar_t>(0.00001);
//   const auto icCenterX = static_cast<scalar_t>(0.0);
//   const auto icCenterY = static_cast<scalar_t>(-0.2);

//   using return_type = PublicProblemEigenMixinCpp<impladvdiff2d::EigenApp<mesh_t>>;
//   return return_type(meshObj,
// 		 problemEnum,
// 		 inviscidFluxRecEnum,
// 		 InviscidFluxScheme::Rusanov,
// 		 viscFluxRecEnum,
// 		 icPulseMagnitude, icSpread, diffusion,
// 		 icCenterX, icCenterY,
// 		 std::move(bcFuncs));

// }
// #endif

// #if !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// template<class mesh_t, class BCsFuncL, class BCsFuncF, class BCsFuncR, class BCsFuncB>
// auto create_periodic_burgers_2d_problem_eigen
// (const mesh_t & meshObj,
//  AdvectionDiffusion2d problemEnum,
//  InviscidFluxReconstruction inviscidFluxRecEnum,
//  ViscousFluxReconstruction viscFluxRecEnum,
//  typename mesh_t::scalar_t icPulseMagnitude,
//  typename mesh_t::scalar_t icSpread,
//  typename mesh_t::scalar_t diffusion,
//  typename mesh_t::scalar_t icCenterX,
//  typename mesh_t::scalar_t icCenterY,
//  BCsFuncL && BCsLeft,
//  BCsFuncF && BCsFront,
//  BCsFuncR && BCsRight,
//  BCsFuncB && BCsBack
//  )
// {

//   using BCFunctorsHolderType = impl::CustomBCsHolder<BCsFuncL, BCsFuncF, BCsFuncR, BCsFuncB>;
//   BCFunctorsHolderType bcFuncs(std::forward<BCsFuncL>(BCsLeft),
// 			       std::forward<BCsFuncF>(BCsFront),
// 			       std::forward<BCsFuncR>(BCsRight),
// 			       std::forward<BCsFuncB>(BCsBack));

//   using return_type = PublicProblemEigenMixinCpp<impladvdiff2d::EigenApp<mesh_t>>;
//   return return_type(meshObj,
// 		 problemEnum,
// 		 inviscidFluxRecEnum,
// 		 InviscidFluxScheme::Rusanov,
// 		 viscFluxRecEnum,
// 		 icPulseMagnitude, icSpread, diffusion,
// 		 icCenterX, icCenterY,
// 		 std::move(bcFuncs));
// }
// #endif

} //end namespace pressiodemoapps
#endif
