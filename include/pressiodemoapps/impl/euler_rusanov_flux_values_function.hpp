/*
//@HEADER
// ************************************************************************
//
// euler_rusanov_flux_values_function.hpp
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

#ifndef PRESSIODEMOAPPS_EE_FLUXES_HPP_
#define PRESSIODEMOAPPS_EE_FLUXES_HPP_

namespace pressiodemoapps{ namespace ee{ namespace impl{

template<class T, typename sc_t>
void ee_rusanov_flux_three_dof(T & F,
			   const T & qL,
			   const T & qR,
			   const sc_t gamma)
{
  constexpr auto one  = static_cast<sc_t>(1);
  constexpr auto two  = static_cast<sc_t>(2);
  constexpr auto half = one/two;

  constexpr sc_t es{1.e-30};
  sc_t FL[3];
  sc_t FR[3];

  // left
  const auto rL = qL(0);
  const auto uL = qL(1)/(rL + es);
  const auto pL = (gamma-1)*( qL(2) - half*rL*(uL*uL) );
  const auto HL = ( qL(2) + pL ) / rL;
  FL[0] = rL*uL;
  FL[1] = rL*uL*uL + pL;
  FL[2] = rL*uL*HL;

  // right
  const auto rR = qR(0);
  const auto uR = qR(1)/(rR + es);
  const auto pR = (gamma-1)*( qR(2) - half*rR*(uR*uR) );
  const auto HR = ( qR(2) + pR ) / rR;
  FR[0] = rR*uR;
  FR[1] = rR*uR*uR + pR;
  FR[2] = rR*uR*HR;

  // compute flux
  const auto RT = std::sqrt(rR/(rL));
  const auto u = (uL+RT*uR)/(1.+RT);
  const auto H = ( HL+RT* HR)/(1.+RT);
  const auto a = std::sqrt( (gamma-1.)*(H - half*(u*u)) );
  const auto smax = std::abs(u) + a;

  F(0) = half*( FL[0] + FR[0] + smax*(qL(0) - qR(0)) );
  F(1) = half*( FL[1] + FR[1] + smax*(qL(1) - qR(1)) );
  F(2) = half*( FL[2] + FR[2] + smax*(qL(2) - qR(2)) );
}

template<typename T, class normal_t, class sc_t>
void ee_rusanov_flux_four_dof(T & F,
			  const T & qL,
			  const T & qR,
			  const normal_t & n,
			  const sc_t gamma)
{
  constexpr auto one  = static_cast<sc_t>(1);
  constexpr auto two  = static_cast<sc_t>(2);
  constexpr auto half = one/two;

  constexpr sc_t es{1.e-30};
  sc_t FL[4];
  sc_t FR[4];

  // left
  const auto rL = qL(0);
  const auto uL = qL(1)/(rL + es);
  const auto vL = qL(2)/(rL + es);
  const auto unL = uL*n[0] + vL*n[1];
  const auto pL = (gamma-1)*( qL(3) - half*rL*(uL*uL+vL*vL) );
  const auto HL = ( qL(3) + pL ) / rL;
  FL[0] = rL*unL;
  FL[1] = rL*unL*uL + pL*n[0];
  FL[2] = rL*unL*vL + pL*n[1];
  FL[3] = rL*unL*HL;

  // right
  const auto rR = qR(0);
  const auto uR = qR(1)/(rR + es);
  const auto vR = qR(2)/(rR + es);
  const auto unR = uR*n[0] + vR*n[1];
  const auto pR = (gamma-1)*( qR(3) - half*rR*(uR*uR+vR*vR) );
  const auto HR = ( qR(3) + pR ) / rR;
  FR[0] = rR*unR;
  FR[1] = rR*unR*uR + pR*n[0];
  FR[2] = rR*unR*vR + pR*n[1];
  FR[3] = rR*unR*HR;

  // compute flux
  const auto RT = std::sqrt(rR/(rL));
  const auto u = (uL+RT*uR)/(1.+RT);
  const auto v = (vL+RT*vR)/(1.+RT);
  const auto H = ( HL+RT* HR)/(1.+RT);
  const auto a = std::sqrt( (gamma-1.)*(H - half*(u*u+v*v)) );
  const auto smax = std::sqrt(u*u+v*v) + a;

  F(0) = half*( FL[0] + FR[0] + smax*(qL(0) - qR(0)) );
  F(1) = half*( FL[1] + FR[1] + smax*(qL(1) - qR(1)) );
  F(2) = half*( FL[2] + FR[2] + smax*(qL(2) - qR(2)) );
  F(3) = half*( FL[3] + FR[3] + smax*(qL(3) - qR(3)) );
}

template<typename T, class normal_t, class sc_t>
void ee_rusanov_flux_five_dof(T & F,
			  const T & qL,
			  const T & qR,
			  const normal_t & n,
			  const sc_t gamma)
{
  constexpr auto one  = static_cast<sc_t>(1);
  constexpr auto two  = static_cast<sc_t>(2);
  constexpr auto half = one/two;

  constexpr sc_t es{1.e-30};
  sc_t FL[5];
  sc_t FR[5];

  // left
  const auto rL = qL(0);
  const auto uL = qL(1)/(rL + es);
  const auto vL = qL(2)/(rL + es);
  const auto wL = qL(3)/(rL + es);
  const auto unL = uL*n[0] + vL*n[1] + wL*n[2];
  const auto pL = (gamma-1)*( qL(4) - half*rL*(uL*uL+vL*vL+wL*wL) );
  const auto HL = ( qL(4) + pL ) / rL;
  FL[0] = rL*unL;
  FL[1] = rL*unL*uL + pL*n[0];
  FL[2] = rL*unL*vL + pL*n[1];
  FL[3] = rL*unL*wL + pL*n[2];
  FL[4] = rL*unL*HL;

  // right
  const auto rR = qR(0);
  const auto uR = qR(1)/(rR + es);
  const auto vR = qR(2)/(rR + es);
  const auto wR = qR(3)/(rR + es);
  const auto unR = uR*n[0] + vR*n[1] + wR*n[2];
  const auto pR = (gamma-1)*( qR(4) - half*rR*(uR*uR+vR*vR+wR*wR) );
  const auto HR = ( qR(4) + pR ) / rR;
  FR[0] = rR*unR;
  FR[1] = rR*unR*uR + pR*n[0];
  FR[2] = rR*unR*vR + pR*n[1];
  FR[3] = rR*unR*wR + pR*n[2];
  FR[4] = rR*unR*HR;

  // compute flux
  const auto RT = std::sqrt(rR/(rL));
  const auto u = (uL+RT*uR)/(1.+RT);
  const auto v = (vL+RT*vR)/(1.+RT);
  const auto w = (wL+RT*wR)/(1.+RT);
  const auto H = ( HL+RT* HR)/(1.+RT);
  const auto a = std::sqrt( (gamma-1.)*(H - half*(u*u + v*v + w*w)) );
  const auto smax = std::sqrt(u*u+v*v+w*w) + a;

  F(0) = half*( FL[0] + FR[0] + smax*(qL(0) - qR(0)) );
  F(1) = half*( FL[1] + FR[1] + smax*(qL(1) - qR(1)) );
  F(2) = half*( FL[2] + FR[2] + smax*(qL(2) - qR(2)) );
  F(3) = half*( FL[3] + FR[3] + smax*(qL(3) - qR(3)) );
  F(4) = half*( FL[4] + FR[4] + smax*(qL(4) - qR(4)) );
}

}}}
#endif
