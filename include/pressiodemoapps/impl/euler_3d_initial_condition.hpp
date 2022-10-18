/*
//@HEADER
// ************************************************************************
//
// euler_3d_initial_condition.hpp
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

#ifndef PRESSIODEMOAPPS_EULER3D_IC_HPP_
#define PRESSIODEMOAPPS_EULER3D_IC_HPP_

#include "../mypi.hpp"

namespace pressiodemoapps{
namespace impleuler3d{

template<class state_type, class mesh_t, class scalar_type>
void euler3dsmoothInitialCondition(state_type & state,
				   const mesh_t & meshObj,
				   const scalar_type gamma)
{
  constexpr int numDofPerCell = 5;
  const auto &x= meshObj.viewX();
  const auto &y= meshObj.viewY();
  const auto &z= meshObj.viewZ();
  const auto gammaMinusOne = gamma - 1.;
  const auto gammaMinusOneInv = 1./gammaMinusOne;

  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {

      const auto ind = i*numDofPerCell;
      prim[0] = 1.0 + 0.2*std::sin(M_PI*(x(i)+y(i)+z(i)));
      prim[1] = 1.0;
      prim[2] = 1.0;
      prim[3] = 1.0;
      prim[4] = 1.;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = prim[0]*prim[2];
      state(ind+3) = prim[0]*prim[3];
      state(ind+4) = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
    }
}

template<class state_type, class mesh_t, class scalar_type>
void sedov3dInitialCondition(state_type & state,
			     const mesh_t & meshObj,
			     const scalar_type gamma)
{
  constexpr int numDofPerCell = 5;

  const auto dx  = meshObj.dx();
  const auto dy  = meshObj.dy();
  const auto dz  = meshObj.dz();
  const auto sRad = 3.*std::min(dx, std::min(dy, dz) );

  const auto &x= meshObj.viewX();
  const auto &y= meshObj.viewY();
  const auto &z= meshObj.viewZ();
  const auto gammaMinusOne = gamma - 1.;
  const auto gammaMinusOneInv = 1./gammaMinusOne;

  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto ind = i*numDofPerCell;
      const auto distX = x(i);
      const auto distY = y(i);
      const auto distZ = z(i);
      const auto xsq = distX*distX;
      const auto ysq = distY*distY;
      const auto zsq = distZ*distZ;
      const auto myR = std::sqrt(xsq+ysq+zsq);

      if (myR <= sRad)
	{
	  prim[0] = 1.0;
	  prim[1] = 0.0;
	  prim[2] = 0.0;
	  prim[3] = 0.0;
	  prim[4] = (3.*gammaMinusOne*0.851072)/(4.*M_PI*sRad*sRad*sRad);
	}

      else
	{
	  prim[0] = 1.;
	  prim[1] = 0.0;
	  prim[2] = 0.0;
	  prim[3] = 0.0;
	  prim[4] = 2.5e-5;
	}

      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = prim[0]*prim[2];
      state(ind+3) = prim[0]*prim[3];
      state(ind+4) = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
    }
}

}}//end namespace

#endif
