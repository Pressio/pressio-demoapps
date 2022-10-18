/*
//@HEADER
// ************************************************************************
//
// diffusion_reaction2d.hpp
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

#ifndef PRESSIODEMOAPPS_DIFFUSION_REACTION_HPP_
#define PRESSIODEMOAPPS_DIFFUSION_REACTION_HPP_

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
enum class DiffusionReaction2d
{
  ProblemA,
  /*
    ds/dt = D (d^2s/d^2x + d^2s/d^2y) + k*s^2 + u(x, y, t)

    BC: ghost cells are set such that s is zero at boundary

    D, k, u(x, y, t) can be provided to constructor
    u(x, y, t) must be a functor:
      void operator()(const scalar_t & x,
		      const scalar_t & y,
		      const scalar_t & time,
		      scalar_t & value);

    Default:
    D, k = 0.01, 0.01
    u(x, y, t) = std::sin(M_PI*x*(y-0.2)) * 4.*std::sin(4.*M_PI*y*x);

   */

  GrayScott
  /*
    two species:
    da/dt = Da Laplacian(a) + F(1-a) - a*b^2
    db/dt = Db Laplacian(v) + (F+k)b + a*b^2

    F = feeding rate, k = kill rate
    Da, Db = diffusion coeffs
    periodic BC
   */
};

// ----------------------------------------------------------
// default source terms
// ----------------------------------------------------------
namespace impldiffusionreaction2d{
template<class scalar_type>
struct DefaultSourceProblemA
{
  void operator()(const scalar_type & x,
		  const scalar_type & y,
		  const scalar_type & evaltime,
		  scalar_type & value)
  {
    (void) evaltime;
    value = std::sin(M_PI*x*(y-0.2)) * 4.*std::sin(4.*M_PI*y*x);
  }
};
}//end namespace impldiffusionreaction1d
}//end namespace pressiodemoapps


// this include is here because needs visiblity of the enums above
#include "./impl/diffusion_reaction_2d_prob_class.hpp"


namespace pressiodemoapps{

// ----------------------------------------------------------
// create default problem
// ----------------------------------------------------------
template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impldiffusionreaction2d::EigenApp<mesh_t>>
  >
RetType
// bindings need unique naming or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
create_diffreac2d_problem_default_for_py
#else
create_problem_eigen
#endif
(const mesh_t & meshObj,
 DiffusionReaction2d problemEnum,
 ViscousFluxReconstruction viscRecEn)
{

  using scalar_t = typename mesh_t::scalar_t;
  if (problemEnum == DiffusionReaction2d::ProblemA)
  {
    return RetType(impldiffusionreaction2d::TagProblemA{},
		   meshObj, viscRecEn,
		   impldiffusionreaction2d::DefaultSourceProblemA<scalar_t>(),
		   0.01, 0.01);
  }

  else if (problemEnum == DiffusionReaction2d::GrayScott)
  {
    return RetType(impldiffusionreaction2d::TagProblemGrayScott{},
		   meshObj, viscRecEn, 0.0002, 0.00005, 0.042, 0.062);
  }

  else{
    throw std::runtime_error("2D diffusion-reaction: invalid problem enum");
  }
}

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impldiffusionreaction2d::EigenApp<mesh_t>>
  >
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
RetType create_diffreac2d_problem_ov1_for_py(const mesh_t & meshObj,
					     DiffusionReaction2d problemEnum)
{
  return create_diffreac2d_problem_default_for_py(meshObj, problemEnum,
						  ViscousFluxReconstruction::FirstOrder);
}
#else
RetType create_problem_eigen(const mesh_t & meshObj,
			     DiffusionReaction2d problemEnum)
{
  return create_problem_eigen(meshObj, problemEnum,
			      ViscousFluxReconstruction::FirstOrder);
}
#endif

// ----------------------------------------------------------
// create problem A with custom diffusion/reaction coeffs
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impldiffusionreaction2d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// bindings need unique naming or we get error associated with overloads
  create_diffreac2d_problem_A_ov1_for_py
#else
  create_diffusion_reaction_2d_problem_A_eigen
#endif
(const mesh_t & meshObj,
 ViscousFluxReconstruction viscFluxRecEnum,
 typename mesh_t::scalar_t diffusion,
 typename mesh_t::scalar_t reaction)
{

  using scalar_t = typename mesh_t::scalar_t;
  return RetType(impldiffusionreaction2d::TagProblemA{},
		 meshObj,
		 viscFluxRecEnum,
		 impldiffusionreaction2d::DefaultSourceProblemA<scalar_t>(),
		 diffusion, reaction);
}

// ----------------------------------------------------------
// create problem A with custom diffusion/reaction and source
// ----------------------------------------------------------
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impldiffusionreaction2d::EigenApp<mesh_t>>
  >
RetType create_diffreac2d_problem_A_ov2_for_py(const mesh_t & meshObj,
					       ViscousFluxReconstruction viscFluxRecEnum,
					       // source term is passed as a Python functor
					       pybind11::object pyFunctor,
					       typename mesh_t::scalar_t diffusion,
					       typename mesh_t::scalar_t reaction)
{
  using scalar_t = typename mesh_t::scalar_t;
  auto sourceWrapper = [=](const scalar_t & x,
			   const scalar_t & y,
			   const scalar_t & evaltime,
			   scalar_t & value){
    value = pyFunctor.attr("__call__")(x, y, evaltime).template cast<scalar_t>();
  };

  return RetType(impldiffusionreaction2d::TagProblemA{},
		 meshObj, viscFluxRecEnum,
		 sourceWrapper, diffusion, reaction);
}

#else

template<
  class mesh_t,
  class SourceType,
  class RetType = PublicProblemEigenMixinCpp<impldiffusionreaction2d::EigenApp<mesh_t>>
  >
RetType create_diffusion_reaction_2d_problem_A_eigen(const mesh_t & meshObj,
						     ViscousFluxReconstruction viscFluxRecEnum,
						     SourceType sourceF,
						     typename mesh_t::scalar_t diffusion,
						     typename mesh_t::scalar_t reaction)
{
  return RetType(impldiffusionreaction2d::TagProblemA{},
		 meshObj, viscFluxRecEnum, sourceF,
		 diffusion, reaction);
}
#endif

// ----------------------------------------------------------
// create gray-scott with custom diffusion/reaction coeffs
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impldiffusionreaction2d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// bindings need unique naming or we get error associated with overloads
  create_diffreac2d_problem_grayscott_ov1_for_py
#else
  create_gray_scott_2d_problem_eigen
#endif
(const mesh_t & meshObj,
 ViscousFluxReconstruction viscFluxRecEnum,
 typename mesh_t::scalar_t diffusion_u,
 typename mesh_t::scalar_t diffusion_v,
 typename mesh_t::scalar_t feedRate,
 typename mesh_t::scalar_t killRate)
{
  return RetType(impldiffusionreaction2d::TagProblemGrayScott{},
		 meshObj,
		 viscFluxRecEnum,
		 diffusion_u, diffusion_v, feedRate, killRate);
}

} //end namespace pressiodemoapps
#endif
