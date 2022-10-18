/*
//@HEADER
// ************************************************************************
//
// euler_compute_energy.hpp
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

#ifndef PRESSIODEMOAPPS_EE_COMPUTE_ENERGY_HPP_
#define PRESSIODEMOAPPS_EE_COMPUTE_ENERGY_HPP_

namespace pressiodemoapps{

/* 1d euler, 3 conserved variables */
template<class scalar_type>
scalar_type eulerEquationsComputeEnergyFromPrimitive2(const scalar_type & gammaMinusOneInv,
						      const std::array<scalar_type, 3> & prim)
{
  const auto one = static_cast<scalar_type>(1);
  const auto two = static_cast<scalar_type>(2);
  const auto half = one/two;

  const auto usq = prim[1]*prim[1];
  return prim[2]*gammaMinusOneInv + half*prim[0]*(usq);
}

template<class scalar_type>
scalar_type eulerEquationsComputeEnergyFromPrimitive(const scalar_type & gamma,
						     const std::array<scalar_type, 3> & prim)
{
  const auto one = static_cast<scalar_type>(1);

  const auto gammaMinusOneInv = one/(gamma - one);
  return eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
}

/* 2d euler, 4 conserved variables */
template<class scalar_type>
scalar_type eulerEquationsComputeEnergyFromPrimitive2(const scalar_type & gammaMinusOneInv,
						      const std::array<scalar_type, 4> & prim)
{
  const auto one = static_cast<scalar_type>(1);
  const auto two = static_cast<scalar_type>(2);
  const auto half = one/two;

  const auto usq = prim[1]*prim[1];
  const auto vsq = prim[2]*prim[2];
  return prim[3]*gammaMinusOneInv + half*prim[0]*(usq+vsq);
}

template<class scalar_type>
scalar_type eulerEquationsComputeEnergyFromPrimitive(const scalar_type & gamma,
						     const std::array<scalar_type, 4> & prim)
{
  const auto one = static_cast<scalar_type>(1);

  const auto gammaMinusOneInv = one/(gamma - one);
  return eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
}

/* 3d euler, 5 conserved variables */
template<class scalar_type>
scalar_type eulerEquationsComputeEnergyFromPrimitive2(const scalar_type & gammaMinusOneInv,
						      const std::array<scalar_type, 5> & prim)
{
  const auto one = static_cast<scalar_type>(1);
  const auto two = static_cast<scalar_type>(2);
  const auto half = one/two;

  const auto usq = prim[1]*prim[1];
  const auto vsq = prim[2]*prim[2];
  const auto wsq = prim[3]*prim[3];
  return prim[4]*gammaMinusOneInv + half*prim[0]*(usq+vsq+wsq);
}

template<class scalar_type>
scalar_type eulerEquationsComputeEnergyFromPrimitive(const scalar_type & gamma,
						     const std::array<scalar_type, 5> & prim)
{
  const auto one = static_cast<scalar_type>(1);

  const auto gammaMinusOneInv = one/(gamma - one);
  return eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
}

}
#endif
