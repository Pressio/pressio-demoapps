/*
//@HEADER
// ************************************************************************
//
// swe_2d_initial_condition.hpp
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

#ifndef PRESSIODEMOAPPS_SWE2D_IC_HPP_
#define PRESSIODEMOAPPS_SWE2D_IC_HPP_

namespace pressiodemoapps{
namespace implswe2d{

template<class state_type, class mesh_t, class scalar_type>
void GaussianPulse(state_type & state,
		   const mesh_t & meshObj,
		   const scalar_type pulseMag,
		   const scalar_type pulseX,
		   const scalar_type pulseY)
{
  constexpr int numDofPerCell = 3;
  const auto &x= meshObj.viewX();
  const auto &y= meshObj.viewY();

  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i){
    const auto ind = i*numDofPerCell;
    const scalar_type dx1 = x(i) - pulseX;
    const scalar_type dy1 = y(i) - pulseY;
    const auto r = std::sqrt(dx1*dx1 + dy1*dy1);
    state(ind)   = static_cast<scalar_type>(1) + pulseMag*std::exp( -(r*r) );
    state(ind+1) = static_cast<scalar_type>(0);
    state(ind+2) = static_cast<scalar_type>(0);
  }
}

template<class state_type, class mesh_t, class scalar_type>
void DoubleGaussianPulse(state_type & state,
			 const mesh_t & meshObj,
			 const scalar_type pulseMag1,
			 const scalar_type pulseX1,
			 const scalar_type pulseY1,
			 const scalar_type pulseMag2,
			 const scalar_type pulseX2,
			 const scalar_type pulseY2)
{
  constexpr int numDofPerCell = 3;
  const auto &x= meshObj.viewX();
  const auto &y= meshObj.viewY();
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
  {
    const auto ind = i*numDofPerCell;
    const scalar_type dx1 = x(i) - pulseX1;
    const scalar_type dy1 = y(i) - pulseY1;
    const auto r1 = std::sqrt(dx1*dx1 + dy1*dy1);

    const scalar_type dx2 = x(i) - pulseX2;
    const scalar_type dy2 = y(i) - pulseY2;
    const auto r2 = std::sqrt(dx2*dx2 + dy2*dy2);
    state(ind) = static_cast<scalar_type>(1) +
		 pulseMag1*std::exp(-(r1*r1)) +
		 pulseMag2*std::exp(-(r2*r2));
    state(ind+1) = static_cast<scalar_type>(0);
    state(ind+2) = static_cast<scalar_type>(0);
  }
}

}}//end namespace
#endif
