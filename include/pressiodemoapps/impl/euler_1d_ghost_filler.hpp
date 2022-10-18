/*
//@HEADER
// ************************************************************************
//
// euler_1d_ghost_filler.hpp
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

#ifndef PRESSIODEMOAPPS_GHOST_FILLER_EE1D_HPP_
#define PRESSIODEMOAPPS_GHOST_FILLER_EE1D_HPP_

namespace pressiodemoapps{ namespace impl{

template<class state_t, class mesh_t, class ghost_t>
class Ghost1dNeumannFiller
{

public:
  Ghost1dNeumannFiller() = delete;
  Ghost1dNeumannFiller(const int stencilSize,
		       const state_t & stateIn,
		       const mesh_t & meshIn,
		       ghost_t & ghostLeft,
		       ghost_t & ghostRight)
    : m_stencilSize(stencilSize),
      m_state(stateIn),
      m_meshObj(meshIn),
      m_ghostLeft(ghostLeft),
      m_ghostRight(ghostRight)
  {}

  template<class index_t>
  void operator()(index_t smPt, int gRow, int ndpc)
  {
    if (ndpc==3){
      threeDofImpl(smPt, gRow);
    }
  }

private:
  template<class index_t>
  void threeDofImpl(index_t smPt, int gRow)
  {
    constexpr int ndpc = 3;
    const auto & graph = m_meshObj.graph();
    const auto cellGID = graph(smPt, 0);
    const auto uIndex  = cellGID*ndpc;

    const auto left0  = graph(smPt, 1);
    const auto right0 = graph(smPt, 2);

    if (left0 == -1){
      m_ghostLeft(gRow, 0) = m_state(uIndex);
      m_ghostLeft(gRow, 1) = m_state(uIndex+1);
      m_ghostLeft(gRow, 2) = m_state(uIndex+2);
    }

    if (right0 == -1){
      m_ghostRight(gRow, 0) = m_state(uIndex);
      m_ghostRight(gRow, 1) = m_state(uIndex+1);
      m_ghostRight(gRow, 2) = m_state(uIndex+2);
    }

    if (m_stencilSize>=5)
      {
	const auto left1  = graph(smPt, 3);
	const auto right1 = graph(smPt, 4);

	if (left1 == -1){
	  const auto ind = right0*ndpc;
	  m_ghostLeft(gRow, 3) = m_state(ind);
	  m_ghostLeft(gRow, 4) = m_state(ind+1);
	  m_ghostLeft(gRow, 5) = m_state(ind+2);
	}

	if (right1 == -1){
	  const auto ind = left0*ndpc;
	  m_ghostRight(gRow, 3) = m_state(ind);
	  m_ghostRight(gRow, 4) = m_state(ind+1);
	  m_ghostRight(gRow, 5) = m_state(ind+2);
	}
      }

    if (m_stencilSize == 7){
      const auto left1  = graph(smPt, 3);
      const auto right1 = graph(smPt, 4);
      const auto left2  = graph(smPt, 5);
      const auto right2 = graph(smPt, 6);

      if (left2 == -1){
	const auto ind = right1*ndpc;
	m_ghostLeft(gRow, 6) = m_state(ind);
	m_ghostLeft(gRow, 7) = m_state(ind+1);
	m_ghostLeft(gRow, 8) = m_state(ind+2);
      }

      if (right2 == -1){
	const auto ind = left1*ndpc;
	m_ghostRight(gRow, 6) = m_state(ind);
	m_ghostRight(gRow, 7) = m_state(ind+1);
	m_ghostRight(gRow, 8) = m_state(ind+2);
      }
    }
  }

private:
  const int m_stencilSize;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  ghost_t & m_ghostLeft;
  ghost_t & m_ghostRight;
};

}}
#endif
