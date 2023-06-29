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

#ifndef PRESSIODEMOAPPS_ADVECTION_DIFFUSION_REACTION_2D_HPP_
#define PRESSIODEMOAPPS_ADVECTION_DIFFUSION_REACTION_2D_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./adapter_cpp.hpp"
#include "./adapter_py.hpp"

namespace pressiodemoapps{

enum class AdvectionDiffusionReaction2d{
  ProblemA
  /*
    d \phi/dt + u dphi/dx + v dphi/dy = D (d^2s/d^2x + d^2s/d^2y) - sigma*\phi + f(x, y, t)

    homog Dirichlet BC

    D, sigma, f(x, y, t) can be passed to constructor

    f(x, y, t) must be a functor:
      void operator()(const scalar_t & x,
		      const scalar_t & y,
		      const scalar_t & time,
		      scalar_t & value);

    Default:
      u_x = 0.5*std::cos( M_PI/3 );
      u_y = 0.5*std::sin( M_PI/3 );
      D = 0.001
      sigma = 1.0
      f(x, y, t) = 1.0
   */
};

// ----------------------------------------------------------
// default source terms
// ----------------------------------------------------------
namespace impladvdiffreac2d{
template<class scalar_type>
struct DefaultSourceProblemA
{
  void operator()(const scalar_type & /*x*/,
		  const scalar_type & /*y*/,
		  const scalar_type & /*evaltime*/,
		  scalar_type & value)
  {
    value = 1.0;
  }
};
}//end namespace impladvdiffreac2d
}//end namespace pressiodemoapps

// this include is here because needs visiblity of the enums above
#include "./impl/advection_diffusion_reaction_2d_prob_class.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// create default problem
// ----------------------------------------------------------
template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impladvdiffreac2d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// bindings need unique naming or we get error associated with overloads
create_advecdiffusionreac2d_problem_default_for_py
#else
create_problem_eigen
#endif
(const mesh_t & meshObj,
 AdvectionDiffusionReaction2d problemEnum,
 InviscidFluxReconstruction inviscidFluxRecEnum,
 ViscousFluxReconstruction viscFluxRecEnum)
{

  using scalar_t = typename mesh_t::scalar_t;
  if (problemEnum == ::pressiodemoapps::AdvectionDiffusionReaction2d::ProblemA)
  {
    // default parameters
    const auto ux = 0.5*std::cos( M_PI/3 );
    const auto uy = 0.5*std::sin( M_PI/3 );
    const auto diffusion = static_cast<scalar_t>(0.001);
    const auto f_source  = static_cast<scalar_t>(1.0);
    const auto sigma_reaction = static_cast<scalar_t>(1.0);

    return RetType(impladvdiffreac2d::TagProblemA{},
		   meshObj,
		   inviscidFluxRecEnum,
		   InviscidFluxScheme::Rusanov,
		   viscFluxRecEnum,
		   ux, uy, diffusion, sigma_reaction,
		   impladvdiffreac2d::DefaultSourceProblemA<scalar_t>());
  }
  else{
    throw std::runtime_error("advection-diffusion2d: invalid problem enum");
  }
}

// ----------------------------------------------------------
// create problem A with custom diffusion, reaction and advection
// ----------------------------------------------------------
template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impladvdiffreac2d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// bindings need unique naming or we get error associated with overloads
create_advecdiffusionreac2d_problem_ov1_for_py
#else
create_problem_eigen
#endif
(const mesh_t & meshObj,
 AdvectionDiffusionReaction2d problemEnum,
 InviscidFluxReconstruction inviscidFluxRecEnum,
 ViscousFluxReconstruction viscFluxRecEnum,
 typename mesh_t::scalar_t ux,
 typename mesh_t::scalar_t uy,
 typename mesh_t::scalar_t diffusion,
 typename mesh_t::scalar_t sigmaReaction)
{

  using scalar_t = typename mesh_t::scalar_t;
  if (problemEnum == ::pressiodemoapps::AdvectionDiffusionReaction2d::ProblemA)
  {
    const auto f_source  = static_cast<scalar_t>(1.0);
    return RetType(impladvdiffreac2d::TagProblemA{},
		   meshObj,
		   inviscidFluxRecEnum,
		   InviscidFluxScheme::Rusanov,
		   viscFluxRecEnum,
		   ux, uy, diffusion, sigmaReaction,
		   impladvdiffreac2d::DefaultSourceProblemA<scalar_t>());
  }
  else{
    throw std::runtime_error("advection-diffusion2d: invalid problem enum");
  }
}

} //end namespace pressiodemoapps
#endif
