/*
//@HEADER
// ************************************************************************
//
// mixin_directional_flux_balance.hpp
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

#ifndef PRESSIODEMOAPPS_MIXIN_DIRECTIONAL_FLUX_BALANCE_NON_TEMP_HPP_
#define PRESSIODEMOAPPS_MIXIN_DIRECTIONAL_FLUX_BALANCE_NON_TEMP_HPP_

namespace pressiodemoapps{ namespace impl{

template<class Parent, class ContainerType, class ScalarType>
struct ComputeDirectionalFluxBalance : Parent
{
private:
  ContainerType & m_V;
  ScalarType m_hInverse;

public:
  template<class ...Args>
  ComputeDirectionalFluxBalance(ContainerType & V,
					   ScalarType hInverse,
					   Args && ...args)
    : Parent(std::forward<Args>(args)...), m_V(V), m_hInverse(hInverse){}

  template<class index_t, class ...Args2>
  void operator()(index_t smPt, int ndpc, Args2 && ...args2)
  {
    Parent::operator()(smPt, ndpc, std::forward<Args2>(args2)...);

    const auto & fluxL = Parent::fluxLeft();
    const auto & fluxR = Parent::fluxRight();
    if (ndpc == 1){
      m_V(smPt) += m_hInverse*(fluxL(0) - fluxR(0));
    }
    else{
      const auto vIndexCurrentCellFirstDof = smPt*ndpc;
      for (int dof=0; dof<ndpc; ++dof){
	m_V(vIndexCurrentCellFirstDof+dof) += m_hInverse*(fluxL(dof) - fluxR(dof));
      }
    }
  }
};

}} //end namespaces
#endif
