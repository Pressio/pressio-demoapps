/*
//@HEADER
// ************************************************************************
//
// functor_reconstruct_from_stencil.hpp
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

#ifndef PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STENCIL_NONTEMP_HPP_
#define PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STENCIL_NONTEMP_HPP_

namespace pressiodemoapps{ namespace impl{

namespace
{
template<class ReconstructedValueType, class StencilDataType>
struct _ReconstructorFromStencilMembers
{
  _ReconstructorFromStencilMembers() = delete;
  _ReconstructorFromStencilMembers(ReconstructionScheme recEn,
				   const StencilDataType & stencilVals,
				   ReconstructedValueType & uMinusHalfNeg,
				   ReconstructedValueType & uMinusHalfPos,
				   ReconstructedValueType & uPlusHalfNeg,
				   ReconstructedValueType & uPlusHalfPos)
    : m_recEn(recEn),
      m_stencilVals(stencilVals),
      m_uMinusHalfNeg(uMinusHalfNeg),
      m_uMinusHalfPos(uMinusHalfPos),
      m_uPlusHalfNeg(uPlusHalfNeg),
      m_uPlusHalfPos(uPlusHalfPos)
  {}

  ::pressiodemoapps::ReconstructionScheme m_recEn;
  const StencilDataType & m_stencilVals;
  ReconstructedValueType & m_uMinusHalfNeg;
  ReconstructedValueType & m_uMinusHalfPos;
  ReconstructedValueType & m_uPlusHalfNeg;
  ReconstructedValueType & m_uPlusHalfPos;
};
} // anonym namespace

template<class ReconstructedValueType, class StencilDataType>
class ReconstructorFromStencil
{
  using members_t = _ReconstructorFromStencilMembers<ReconstructedValueType, StencilDataType>;
  members_t m_memb;

public:
  template<typename ...Args>
  ReconstructorFromStencil(Args && ... args)
    : m_memb(std::forward<Args>(args)...){}

  ReconstructionScheme reconstructionScheme() const{
    return m_memb.m_recEn;
  }
  const ReconstructedValueType & reconstructionLeftNeg() const{
    return m_memb.m_uMinusHalfNeg;
  }
  const ReconstructedValueType & reconstructionLeftPos() const{
    return m_memb.m_uMinusHalfPos;
  }
  const ReconstructedValueType & reconstructionRightNeg() const{
    return m_memb.m_uPlusHalfNeg;
  }
  const ReconstructedValueType & reconstructionRightPos() const{
    return m_memb.m_uPlusHalfPos;
  }

  template<class IndexType>
  void operator()(IndexType /*unused*/, int ndpc)
  {
    if (ndpc == 1){
      oneDofPerCellImpl();
    }

    else if (ndpc == 2){
      twoDofPerCellImpl();
    }

    else if (ndpc == 3){
      threeDofPerCellImpl();
    }

    else if (ndpc == 4){
      fourDofPerCellImpl();
    }

    else if (ndpc == 5){
      fiveDofPerCellImpl();
    }
  }

private:
  void oneDofPerCellImpl()
  {
    switch(m_memb.m_recEn){
    case ::pressiodemoapps::ReconstructionScheme::FirstOrder:
      {
	m_memb.m_uMinusHalfNeg(0) = m_memb.m_stencilVals(0);
	m_memb.m_uMinusHalfPos(0) = m_memb.m_stencilVals(1);
	m_memb.m_uPlusHalfNeg(0)  = m_memb.m_stencilVals(1);
	m_memb.m_uPlusHalfPos(0)  = m_memb.m_stencilVals(2);
	break;
      }

    case ::pressiodemoapps::ReconstructionScheme::Weno3:
      {
	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(0),
			       m_memb.m_uMinusHalfPos(0),
			       m_memb.m_uPlusHalfNeg(0),
			       m_memb.m_uPlusHalfPos(0),
			       m_memb.m_stencilVals);
	break;
      }

    case ::pressiodemoapps::ReconstructionScheme::Weno5:
      {
	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(0),
			       m_memb.m_uMinusHalfPos(0),
			       m_memb.m_uPlusHalfNeg(0),
			       m_memb.m_uPlusHalfPos(0),
			       m_memb.m_stencilVals);
	break;
      }
    }
  }

  void twoDofPerCellImpl()
  {
    switch(m_memb.m_recEn)
      {
      case ::pressiodemoapps::ReconstructionScheme::FirstOrder:{
	m_memb.m_uMinusHalfNeg(0) = m_memb.m_stencilVals(0);
	m_memb.m_uMinusHalfPos(0) = m_memb.m_stencilVals(2);
	m_memb.m_uPlusHalfNeg(0)  = m_memb.m_stencilVals(2);
	m_memb.m_uPlusHalfPos(0)  = m_memb.m_stencilVals(4);

	m_memb.m_uMinusHalfNeg(1) = m_memb.m_stencilVals(1);
	m_memb.m_uMinusHalfPos(1) = m_memb.m_stencilVals(3);
	m_memb.m_uPlusHalfNeg(1)  = m_memb.m_stencilVals(3);
	m_memb.m_uPlusHalfPos(1)  = m_memb.m_stencilVals(5);
	break;
      }

      case ::pressiodemoapps::ReconstructionScheme::Weno3:{
	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(0),
			       m_memb.m_uMinusHalfPos(0),
			       m_memb.m_uPlusHalfNeg(0),
			       m_memb.m_uPlusHalfPos(0),
			       m_memb.m_stencilVals(0),
			       m_memb.m_stencilVals(2),
			       m_memb.m_stencilVals(4),
			       m_memb.m_stencilVals(6),
			       m_memb.m_stencilVals(8));

	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(1),
			       m_memb.m_uMinusHalfPos(1),
			       m_memb.m_uPlusHalfNeg(1),
			       m_memb.m_uPlusHalfPos(1),
			       m_memb.m_stencilVals(1),
			       m_memb.m_stencilVals(3),
			       m_memb.m_stencilVals(5),
			       m_memb.m_stencilVals(7),
			       m_memb.m_stencilVals(9));
	break;
      }

      case ::pressiodemoapps::ReconstructionScheme::Weno5:{
	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(0),
			       m_memb.m_uMinusHalfPos(0),
			       m_memb.m_uPlusHalfNeg(0),
			       m_memb.m_uPlusHalfPos(0),
			       m_memb.m_stencilVals(0),
			       m_memb.m_stencilVals(2),
			       m_memb.m_stencilVals(4),
			       m_memb.m_stencilVals(6),
			       m_memb.m_stencilVals(8),
			       m_memb.m_stencilVals(10),
			       m_memb.m_stencilVals(12));

	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(1),
			       m_memb.m_uMinusHalfPos(1),
			       m_memb.m_uPlusHalfNeg(1),
			       m_memb.m_uPlusHalfPos(1),
			       m_memb.m_stencilVals(1),
			       m_memb.m_stencilVals(3),
			       m_memb.m_stencilVals(5),
			       m_memb.m_stencilVals(7),
			       m_memb.m_stencilVals(9),
			       m_memb.m_stencilVals(11),
			       m_memb.m_stencilVals(13));
	break;
      }
      }
  }

  void threeDofPerCellImpl()
  {
    switch(m_memb.m_recEn)
      {
      case ::pressiodemoapps::ReconstructionScheme::FirstOrder:{
	m_memb.m_uMinusHalfNeg(0) = m_memb.m_stencilVals(0);
	m_memb.m_uMinusHalfPos(0) = m_memb.m_stencilVals(3);
	m_memb.m_uPlusHalfNeg(0)  = m_memb.m_stencilVals(3);
	m_memb.m_uPlusHalfPos(0)  = m_memb.m_stencilVals(6);

	m_memb.m_uMinusHalfNeg(1) = m_memb.m_stencilVals(1);
	m_memb.m_uMinusHalfPos(1) = m_memb.m_stencilVals(4);
	m_memb.m_uPlusHalfNeg(1)  = m_memb.m_stencilVals(4);
	m_memb.m_uPlusHalfPos(1)  = m_memb.m_stencilVals(7);

	m_memb.m_uMinusHalfNeg(2) = m_memb.m_stencilVals(2);
	m_memb.m_uMinusHalfPos(2) = m_memb.m_stencilVals(5);
	m_memb.m_uPlusHalfNeg(2)  = m_memb.m_stencilVals(5);
	m_memb.m_uPlusHalfPos(2)  = m_memb.m_stencilVals(8);
	break;
      }

      case ::pressiodemoapps::ReconstructionScheme::Weno3:{
	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(0),
			       m_memb.m_uMinusHalfPos(0),
			       m_memb.m_uPlusHalfNeg(0),
			       m_memb.m_uPlusHalfPos(0),
			       m_memb.m_stencilVals(0),
			       m_memb.m_stencilVals(3),
			       m_memb.m_stencilVals(6),
			       m_memb.m_stencilVals(9),
			       m_memb.m_stencilVals(12));

	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(1),
			       m_memb.m_uMinusHalfPos(1),
			       m_memb.m_uPlusHalfNeg(1),
			       m_memb.m_uPlusHalfPos(1),
			       m_memb.m_stencilVals(1),
			       m_memb.m_stencilVals(4),
			       m_memb.m_stencilVals(7),
			       m_memb.m_stencilVals(10),
			       m_memb.m_stencilVals(13));

	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(2),
			       m_memb.m_uMinusHalfPos(2),
			       m_memb.m_uPlusHalfNeg(2),
			       m_memb.m_uPlusHalfPos(2),
			       m_memb.m_stencilVals(2),
			       m_memb.m_stencilVals(5),
			       m_memb.m_stencilVals(8),
			       m_memb.m_stencilVals(11),
			       m_memb.m_stencilVals(14));
	break;
      }

      case ::pressiodemoapps::ReconstructionScheme::Weno5:{
	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(0),
			       m_memb.m_uMinusHalfPos(0),
			       m_memb.m_uPlusHalfNeg(0),
			       m_memb.m_uPlusHalfPos(0),
			       m_memb.m_stencilVals(0),
			       m_memb.m_stencilVals(3),
			       m_memb.m_stencilVals(6),
			       m_memb.m_stencilVals(9),
			       m_memb.m_stencilVals(12),
			       m_memb.m_stencilVals(15),
			       m_memb.m_stencilVals(18));

	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(1),
			       m_memb.m_uMinusHalfPos(1),
			       m_memb.m_uPlusHalfNeg(1),
			       m_memb.m_uPlusHalfPos(1),
			       m_memb.m_stencilVals(1),
			       m_memb.m_stencilVals(4),
			       m_memb.m_stencilVals(7),
			       m_memb.m_stencilVals(10),
			       m_memb.m_stencilVals(13),
			       m_memb.m_stencilVals(16),
			       m_memb.m_stencilVals(19));

	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(2),
			       m_memb.m_uMinusHalfPos(2),
			       m_memb.m_uPlusHalfNeg(2),
			       m_memb.m_uPlusHalfPos(2),
			       m_memb.m_stencilVals(2),
			       m_memb.m_stencilVals(5),
			       m_memb.m_stencilVals(8),
			       m_memb.m_stencilVals(11),
			       m_memb.m_stencilVals(14),
			       m_memb.m_stencilVals(17),
			       m_memb.m_stencilVals(20));
	break;
      }
      }
  }

  void fourDofPerCellImpl()
  {
    switch(m_memb.m_recEn)
      {
      case ::pressiodemoapps::ReconstructionScheme::FirstOrder:{
	m_memb.m_uMinusHalfNeg(0) = m_memb.m_stencilVals(0);
	m_memb.m_uMinusHalfPos(0) = m_memb.m_stencilVals(4);
	m_memb.m_uPlusHalfNeg(0)  = m_memb.m_stencilVals(4);
	m_memb.m_uPlusHalfPos(0)  = m_memb.m_stencilVals(8);

	m_memb.m_uMinusHalfNeg(1) = m_memb.m_stencilVals(1);
	m_memb.m_uMinusHalfPos(1) = m_memb.m_stencilVals(5);
	m_memb.m_uPlusHalfNeg(1)  = m_memb.m_stencilVals(5);
	m_memb.m_uPlusHalfPos(1)  = m_memb.m_stencilVals(9);

	m_memb.m_uMinusHalfNeg(2) = m_memb.m_stencilVals(2);
	m_memb.m_uMinusHalfPos(2) = m_memb.m_stencilVals(6);
	m_memb.m_uPlusHalfNeg(2)  = m_memb.m_stencilVals(6);
	m_memb.m_uPlusHalfPos(2)  = m_memb.m_stencilVals(10);

	m_memb.m_uMinusHalfNeg(3) = m_memb.m_stencilVals(3);
	m_memb.m_uMinusHalfPos(3) = m_memb.m_stencilVals(7);
	m_memb.m_uPlusHalfNeg(3)  = m_memb.m_stencilVals(7);
	m_memb.m_uPlusHalfPos(3)  = m_memb.m_stencilVals(11);
	break;
      }

      case ::pressiodemoapps::ReconstructionScheme::Weno3:{
	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(0),
			       m_memb.m_uMinusHalfPos(0),
			       m_memb.m_uPlusHalfNeg(0),
			       m_memb.m_uPlusHalfPos(0),
			       m_memb.m_stencilVals(0),
			       m_memb.m_stencilVals(4),
			       m_memb.m_stencilVals(8),
			       m_memb.m_stencilVals(12),
			       m_memb.m_stencilVals(16));

	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(1),
			       m_memb.m_uMinusHalfPos(1),
			       m_memb.m_uPlusHalfNeg(1),
			       m_memb.m_uPlusHalfPos(1),
			       m_memb.m_stencilVals(1),
			       m_memb.m_stencilVals(5),
			       m_memb.m_stencilVals(9),
			       m_memb.m_stencilVals(13),
			       m_memb.m_stencilVals(17));

	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(2),
			       m_memb.m_uMinusHalfPos(2),
			       m_memb.m_uPlusHalfNeg(2),
			       m_memb.m_uPlusHalfPos(2),
			       m_memb.m_stencilVals(2),
			       m_memb.m_stencilVals(6),
			       m_memb.m_stencilVals(10),
			       m_memb.m_stencilVals(14),
			       m_memb.m_stencilVals(18));

	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(3),
			       m_memb.m_uMinusHalfPos(3),
			       m_memb.m_uPlusHalfNeg(3),
			       m_memb.m_uPlusHalfPos(3),
			       m_memb.m_stencilVals(3),
			       m_memb.m_stencilVals(7),
			       m_memb.m_stencilVals(11),
			       m_memb.m_stencilVals(15),
			       m_memb.m_stencilVals(19));
	break;
      }

      case ::pressiodemoapps::ReconstructionScheme::Weno5:{
	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(0),
			       m_memb.m_uMinusHalfPos(0),
			       m_memb.m_uPlusHalfNeg(0),
			       m_memb.m_uPlusHalfPos(0),
			       m_memb.m_stencilVals(0),
			       m_memb.m_stencilVals(4),
			       m_memb.m_stencilVals(8),
			       m_memb.m_stencilVals(12),
			       m_memb.m_stencilVals(16),
			       m_memb.m_stencilVals(20),
			       m_memb.m_stencilVals(24));

	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(1),
			       m_memb.m_uMinusHalfPos(1),
			       m_memb.m_uPlusHalfNeg(1),
			       m_memb.m_uPlusHalfPos(1),
			       m_memb.m_stencilVals(1),
			       m_memb.m_stencilVals(5),
			       m_memb.m_stencilVals(9),
			       m_memb.m_stencilVals(13),
			       m_memb.m_stencilVals(17),
			       m_memb.m_stencilVals(21),
			       m_memb.m_stencilVals(25));

	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(2),
			       m_memb.m_uMinusHalfPos(2),
			       m_memb.m_uPlusHalfNeg(2),
			       m_memb.m_uPlusHalfPos(2),
			       m_memb.m_stencilVals(2),
			       m_memb.m_stencilVals(6),
			       m_memb.m_stencilVals(10),
			       m_memb.m_stencilVals(14),
			       m_memb.m_stencilVals(18),
			       m_memb.m_stencilVals(22),
			       m_memb.m_stencilVals(26));

	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(3),
			       m_memb.m_uMinusHalfPos(3),
			       m_memb.m_uPlusHalfNeg(3),
			       m_memb.m_uPlusHalfPos(3),
			       m_memb.m_stencilVals(3),
			       m_memb.m_stencilVals(7),
			       m_memb.m_stencilVals(11),
			       m_memb.m_stencilVals(15),
			       m_memb.m_stencilVals(19),
			       m_memb.m_stencilVals(23),
			       m_memb.m_stencilVals(27));
	break;
      }
      }
  }

  void fiveDofPerCellImpl()
  {
    switch(m_memb.m_recEn){
    case ::pressiodemoapps::ReconstructionScheme::FirstOrder:{
      m_memb.m_uMinusHalfNeg(0) = m_memb.m_stencilVals(0);
      m_memb.m_uMinusHalfPos(0) = m_memb.m_stencilVals(5);
      m_memb.m_uPlusHalfNeg(0)  = m_memb.m_stencilVals(5);
      m_memb.m_uPlusHalfPos(0)  = m_memb.m_stencilVals(10);

      m_memb.m_uMinusHalfNeg(1) = m_memb.m_stencilVals(1);
      m_memb.m_uMinusHalfPos(1) = m_memb.m_stencilVals(6);
      m_memb.m_uPlusHalfNeg(1)  = m_memb.m_stencilVals(6);
      m_memb.m_uPlusHalfPos(1)  = m_memb.m_stencilVals(11);

      m_memb.m_uMinusHalfNeg(2) = m_memb.m_stencilVals(2);
      m_memb.m_uMinusHalfPos(2) = m_memb.m_stencilVals(7);
      m_memb.m_uPlusHalfNeg(2)  = m_memb.m_stencilVals(7);
      m_memb.m_uPlusHalfPos(2)  = m_memb.m_stencilVals(12);

      m_memb.m_uMinusHalfNeg(3) = m_memb.m_stencilVals(3);
      m_memb.m_uMinusHalfPos(3) = m_memb.m_stencilVals(8);
      m_memb.m_uPlusHalfNeg(3)  = m_memb.m_stencilVals(8);
      m_memb.m_uPlusHalfPos(3)  = m_memb.m_stencilVals(13);

      m_memb.m_uMinusHalfNeg(4) = m_memb.m_stencilVals(4);
      m_memb.m_uMinusHalfPos(4) = m_memb.m_stencilVals(9);
      m_memb.m_uPlusHalfNeg(4)  = m_memb.m_stencilVals(9);
      m_memb.m_uPlusHalfPos(4)  = m_memb.m_stencilVals(14);
      break;
    }

    case ::pressiodemoapps::ReconstructionScheme::Weno3:{
      pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(0),
			     m_memb.m_uMinusHalfPos(0),
			     m_memb.m_uPlusHalfNeg(0),
			     m_memb.m_uPlusHalfPos(0),
			     m_memb.m_stencilVals(0),
			     m_memb.m_stencilVals(5),
			     m_memb.m_stencilVals(10),
			     m_memb.m_stencilVals(15),
			     m_memb.m_stencilVals(20));

      pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(1),
			     m_memb.m_uMinusHalfPos(1),
			     m_memb.m_uPlusHalfNeg(1),
			     m_memb.m_uPlusHalfPos(1),
			     m_memb.m_stencilVals(1),
			     m_memb.m_stencilVals(6),
			     m_memb.m_stencilVals(11),
			     m_memb.m_stencilVals(16),
			     m_memb.m_stencilVals(21));

      pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(2),
			     m_memb.m_uMinusHalfPos(2),
			     m_memb.m_uPlusHalfNeg(2),
			     m_memb.m_uPlusHalfPos(2),
			     m_memb.m_stencilVals(2),
			     m_memb.m_stencilVals(7),
			     m_memb.m_stencilVals(12),
			     m_memb.m_stencilVals(17),
			     m_memb.m_stencilVals(22));

      pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(3),
			     m_memb.m_uMinusHalfPos(3),
			     m_memb.m_uPlusHalfNeg(3),
			     m_memb.m_uPlusHalfPos(3),
			     m_memb.m_stencilVals(3),
			     m_memb.m_stencilVals(8),
			     m_memb.m_stencilVals(13),
			     m_memb.m_stencilVals(18),
			     m_memb.m_stencilVals(23));

      pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(4),
			     m_memb.m_uMinusHalfPos(4),
			     m_memb.m_uPlusHalfNeg(4),
			     m_memb.m_uPlusHalfPos(4),
			     m_memb.m_stencilVals(4),
			     m_memb.m_stencilVals(9),
			     m_memb.m_stencilVals(14),
			     m_memb.m_stencilVals(19),
			     m_memb.m_stencilVals(24));
      break;
    }

    case ::pressiodemoapps::ReconstructionScheme::Weno5:{
      pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(0),
			     m_memb.m_uMinusHalfPos(0),
			     m_memb.m_uPlusHalfNeg(0),
			     m_memb.m_uPlusHalfPos(0),
			     m_memb.m_stencilVals(0),
			     m_memb.m_stencilVals(5),
			     m_memb.m_stencilVals(10),
			     m_memb.m_stencilVals(15),
			     m_memb.m_stencilVals(20),
			     m_memb.m_stencilVals(25),
			     m_memb.m_stencilVals(30));

      pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(1),
			     m_memb.m_uMinusHalfPos(1),
			     m_memb.m_uPlusHalfNeg(1),
			     m_memb.m_uPlusHalfPos(1),
			     m_memb.m_stencilVals(1),
			     m_memb.m_stencilVals(6),
			     m_memb.m_stencilVals(11),
			     m_memb.m_stencilVals(16),
			     m_memb.m_stencilVals(21),
			     m_memb.m_stencilVals(26),
			     m_memb.m_stencilVals(31));

      pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(2),
			     m_memb.m_uMinusHalfPos(2),
			     m_memb.m_uPlusHalfNeg(2),
			     m_memb.m_uPlusHalfPos(2),
			     m_memb.m_stencilVals(2),
			     m_memb.m_stencilVals(7),
			     m_memb.m_stencilVals(12),
			     m_memb.m_stencilVals(17),
			     m_memb.m_stencilVals(22),
			     m_memb.m_stencilVals(27),
			     m_memb.m_stencilVals(32));

      pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(3),
			     m_memb.m_uMinusHalfPos(3),
			     m_memb.m_uPlusHalfNeg(3),
			     m_memb.m_uPlusHalfPos(3),
			     m_memb.m_stencilVals(3),
			     m_memb.m_stencilVals(8),
			     m_memb.m_stencilVals(13),
			     m_memb.m_stencilVals(18),
			     m_memb.m_stencilVals(23),
			     m_memb.m_stencilVals(28),
			     m_memb.m_stencilVals(33));

      pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(4),
			     m_memb.m_uMinusHalfPos(4),
			     m_memb.m_uPlusHalfNeg(4),
			     m_memb.m_uPlusHalfPos(4),
			     m_memb.m_stencilVals(4),
			     m_memb.m_stencilVals(13),
			     m_memb.m_stencilVals(18),
			     m_memb.m_stencilVals(23),
			     m_memb.m_stencilVals(28),
			     m_memb.m_stencilVals(33),
			     m_memb.m_stencilVals(38));
      break;
    }
    }
  }

};

}}
#endif
