/*
//@HEADER
// ************************************************************************
//
// advection_diffusion_2d_flux_mixin.hpp
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

#ifndef PRESSIODEMOAPPS_BURGERS2D_FLUX_FUNCTOR_HPP_
#define PRESSIODEMOAPPS_BURGERS2D_FLUX_FUNCTOR_HPP_

namespace pressiodemoapps{ namespace impladvdiff2d{

template<class Parent, class scalar_type, class flux_type>
struct ComputeDirectionalFluxValues : Parent
{
private:
  const std::array<scalar_type, 2> m_normal;
  ::pressiodemoapps::AdvectionDiffusion2d m_probEn;
  InviscidFluxScheme m_fluxEnum;
  flux_type & m_fluxL;
  flux_type & m_fluxR;

public:
  template<class ...Args>
  ComputeDirectionalFluxValues(::pressiodemoapps::AdvectionDiffusion2d probEn,
			       InviscidFluxScheme fluxEnum,
			       const std::array<scalar_type, 2> normal,
			       flux_type & fluxL,
			       flux_type & fluxR,
			       Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_normal(normal),
      m_probEn(probEn),
      m_fluxEnum(fluxEnum),
      m_fluxL(fluxL),
      m_fluxR(fluxR)
  {}

  const flux_type & fluxLeft()  const { return m_fluxL; }
  const flux_type & fluxRight() const { return m_fluxR; }

  template<class IndexType>
  void operator()(IndexType smPt, int ndpc)
  {
    Parent::operator()(smPt, ndpc);

    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();

    switch(m_fluxEnum)
      {
      case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	burgers_rusanov_flux_2d<scalar_type>(m_fluxL, uMinusHalfNeg,
					     uMinusHalfPos, m_normal);
	burgers_rusanov_flux_2d<scalar_type>(m_fluxR, uPlusHalfNeg,
					     uPlusHalfPos, m_normal);
	break;
      }
  }
};

template<class Parent, class scalar_type, class flux_jac_type>
struct ComputeDirectionalFluxJacobians : Parent
{
private:
  const std::array<scalar_type, 2> m_normal;
  ::pressiodemoapps::AdvectionDiffusion2d m_probEn;
  InviscidFluxScheme m_fluxEnum;
  flux_jac_type & m_fluxJacLNeg;
  flux_jac_type & m_fluxJacLPos;
  flux_jac_type & m_fluxJacRNeg;
  flux_jac_type & m_fluxJacRPos;

public:
  template<class ...Args>
  ComputeDirectionalFluxJacobians(::pressiodemoapps::AdvectionDiffusion2d probEn,
				  InviscidFluxScheme fluxEnum,
				  const std::array<scalar_type, 2> normal,
				  flux_jac_type & fluxJacLNeg,
				  flux_jac_type & fluxJacLPos,
				  flux_jac_type & fluxJacRNeg,
				  flux_jac_type & fluxJacRPos,
				  Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_normal(normal),
      m_probEn(probEn),
      m_fluxEnum(fluxEnum),
      m_fluxJacLNeg(fluxJacLNeg), m_fluxJacLPos(fluxJacLPos),
      m_fluxJacRNeg(fluxJacRNeg), m_fluxJacRPos(fluxJacRPos)
  {}

  const flux_jac_type & fluxJacLNeg()  const { return m_fluxJacLNeg; }
  const flux_jac_type & fluxJacLPos()  const { return m_fluxJacLPos; }
  const flux_jac_type & fluxJacRNeg()  const { return m_fluxJacRNeg; }
  const flux_jac_type & fluxJacRPos()  const { return m_fluxJacRPos; }

  template<class IndexType>
  void operator()(IndexType smPt, int ndpc)
  {
    Parent::operator()(smPt, ndpc);

    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();

    switch(m_fluxEnum)
      {

      case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	burgers_rusanov_flux_jacobian_2d<scalar_type>(m_fluxJacLNeg, m_fluxJacLPos,
						      uMinusHalfNeg, uMinusHalfPos,
						      m_normal);
	burgers_rusanov_flux_jacobian_2d<scalar_type>(m_fluxJacRNeg, m_fluxJacRPos,
						      uPlusHalfNeg, uPlusHalfPos,
						      m_normal);
	break;
      }
  }
};

template<class Parent, class scalar_type, class flux_type, class flux_jac_type>
struct ComputeDirectionalFluxValuesAndJacobians : Parent
{
private:
  const std::array<scalar_type, 2> m_normal;
  ::pressiodemoapps::AdvectionDiffusion2d m_probEn;
  InviscidFluxScheme m_fluxEnum;
  flux_type & m_fluxL;
  flux_type & m_fluxR;
  flux_jac_type & m_fluxJacLNeg;
  flux_jac_type & m_fluxJacLPos;
  flux_jac_type & m_fluxJacRNeg;
  flux_jac_type & m_fluxJacRPos;

public:
  template<class ...Args>
  ComputeDirectionalFluxValuesAndJacobians(::pressiodemoapps::AdvectionDiffusion2d probEn,
					   InviscidFluxScheme fluxEnum,
					   const std::array<scalar_type, 2> normal,
					   flux_type & fluxL,
					   flux_type & fluxR,
					   flux_jac_type & fluxJacLNeg,
					   flux_jac_type & fluxJacLPos,
					   flux_jac_type & fluxJacRNeg,
					   flux_jac_type & fluxJacRPos,
					   Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_normal(normal),
      m_probEn(probEn),
      m_fluxEnum(fluxEnum),
      m_fluxL(fluxL), m_fluxR(fluxR),
      m_fluxJacLNeg(fluxJacLNeg), m_fluxJacLPos(fluxJacLPos),
      m_fluxJacRNeg(fluxJacRNeg), m_fluxJacRPos(fluxJacRPos)
  {}

  const flux_type & fluxLeft()  const { return m_fluxL; }
  const flux_type & fluxRight() const { return m_fluxR; }
  const flux_jac_type & fluxJacLNeg()  const { return m_fluxJacLNeg; }
  const flux_jac_type & fluxJacLPos()  const { return m_fluxJacLPos; }
  const flux_jac_type & fluxJacRNeg()  const { return m_fluxJacRNeg; }
  const flux_jac_type & fluxJacRPos()  const { return m_fluxJacRPos; }

  template<class IndexType>
  void operator()(IndexType smPt, int ndpc)
  {
    Parent::operator()(smPt, ndpc);

    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();

    switch(m_fluxEnum)
      {
      case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	burgers_rusanov_flux_2d<scalar_type>(m_fluxL, uMinusHalfNeg, uMinusHalfPos, m_normal);
	burgers_rusanov_flux_2d<scalar_type>(m_fluxR, uPlusHalfNeg, uPlusHalfPos, m_normal);

	burgers_rusanov_flux_jacobian_2d<scalar_type>(m_fluxJacLNeg, m_fluxJacLPos,
						      uMinusHalfNeg, uMinusHalfPos,
						      m_normal);
	burgers_rusanov_flux_jacobian_2d<scalar_type>(m_fluxJacRNeg, m_fluxJacRPos,
						      uPlusHalfNeg, uPlusHalfPos,
						      m_normal);
	break;
      }
  }
};

}}
#endif
