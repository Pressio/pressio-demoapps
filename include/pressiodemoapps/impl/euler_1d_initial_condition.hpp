/*
//@HEADER
// ************************************************************************
//
// euler_1d_initial_condition.hpp
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

#ifndef PRESSIODEMOAPPS_EE1D_IC_HPP_
#define PRESSIODEMOAPPS_EE1D_IC_HPP_

#include "../mypi.hpp"

namespace pressiodemoapps{ namespace impleuler1d{

template<class state_type, class mesh_t, class scalar_type>
void euler1dsineInitialCondition(state_type & state,
				 const mesh_t & meshObj,
				 const scalar_type gamma)
{
  constexpr int numDofPerCell = 3;
  constexpr auto one  = static_cast<scalar_type>(1);

  const auto & x= meshObj.viewX();
  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      prim[0] = one + static_cast<scalar_type>(0.2)*std::sin(M_PI*x(i));
      prim[1] = one;
      prim[2] = one;

      const auto ind = i*numDofPerCell;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = eulerEquationsComputeEnergyFromPrimitive(gamma, prim);
    }
}

template<class state_type, class mesh_t, class scalar_type>
void sod1dInitialCondition(state_type & state,
			   const mesh_t & meshObj,
			   const scalar_type gamma)
{

  constexpr int numDofPerCell = 3;
  constexpr auto zero = static_cast<scalar_type>(0);
  constexpr auto one  = static_cast<scalar_type>(1);

  const auto & x= meshObj.viewX();
  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      if (x(i) <= zero){
	prim[0] = one;
	prim[1] = zero;
	prim[2] = one;
      }

      if (x(i) > zero){
	prim[0] = static_cast<scalar_type>(0.125);
	prim[1] = zero;
	prim[2] = static_cast<scalar_type>(0.1);
      }

      const auto ind = i*numDofPerCell;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = eulerEquationsComputeEnergyFromPrimitive(gamma, prim);
    }
}

template<class state_type, class mesh_t, class scalar_type>
void lax1dInitialCondition(state_type & state,
			   const mesh_t & meshObj,
			   const scalar_type gamma)
{

  constexpr int numDofPerCell = 3;
  constexpr auto zero = static_cast<scalar_type>(0);

  const auto & x= meshObj.viewX();
  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto ind = i*numDofPerCell;
      if (x(i) <= zero){
	prim[0] = static_cast<scalar_type>(0.445);
	prim[1] = static_cast<scalar_type>(0.698);
	prim[2] = static_cast<scalar_type>(3.528);
      }
      else if (x(i) > zero){
	prim[0] = static_cast<scalar_type>(0.5);
	prim[1] = zero;
	prim[2] = static_cast<scalar_type>(0.571);
      }

      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = eulerEquationsComputeEnergyFromPrimitive(gamma, prim);
    }
}

template<class state_type, class mesh_t, class scalar_type>
void shuOsherInitialCondition(state_type & state,
			   const mesh_t & meshObj,
			   const scalar_type gamma)
{

  constexpr int numDofPerCell = 3;
  constexpr auto zero    = static_cast<scalar_type>(0);
  constexpr auto one     = static_cast<scalar_type>(1);
  constexpr auto three   = static_cast<scalar_type>(3);
  constexpr auto four    = static_cast<scalar_type>(4);
  constexpr auto negFour = -four;
  constexpr auto five    = static_cast<scalar_type>(5);
  constexpr auto seven   = static_cast<scalar_type>(7);
  constexpr auto twentySeven  = static_cast<scalar_type>(27);
  constexpr auto thirtyOne  = static_cast<scalar_type>(31);

  const auto & x= meshObj.viewX();
  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto ind = i*numDofPerCell;
      if (x(i) <= negFour){
	prim[0] = twentySeven/seven;
	prim[1] = static_cast<scalar_type>(2.629369);
	prim[2] = thirtyOne/three;
      }
      else{
	prim[0] = one + (one/five)*std::sin(five*x(i));
	prim[1] = zero;
	prim[2] = one;
      }

      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = eulerEquationsComputeEnergyFromPrimitive(gamma, prim);
    }
}

}}//end namespace
#endif
