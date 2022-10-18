/*
//@HEADER
// ************************************************************************
//
// diffusion_reaction_2d_ghost_filler.hpp
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

#ifndef PRESSIODEMOAPPS_GHOST_FILLER_DIFFREAC_2D_HPP_
#define PRESSIODEMOAPPS_GHOST_FILLER_DIFFREAC_2D_HPP_

namespace pressiodemoapps{
namespace impldiffusionreaction2d{

template<class state_t, class mesh_t, class ghost_t>
class GhostFillerProblemA2d
{

public:
  GhostFillerProblemA2d() = delete;
  GhostFillerProblemA2d(const int stencilSize,
			const state_t & stateIn,
			const mesh_t & meshIn,
			ghost_t & ghostLeft,
			ghost_t & ghostFront,
			ghost_t & ghostRight,
			ghost_t & ghostBack)
    : m_stencilSize(stencilSize),
      m_state(stateIn),
      m_meshObj(meshIn),
      m_ghostLeft(ghostLeft),
      m_ghostFront(ghostFront),
      m_ghostRight(ghostRight),
      m_ghostBack(ghostBack)
  {}

  template<class index_t>
  void operator()(index_t smPt, int gRow)
  {
    const auto & graph = m_meshObj.graph();
    const auto cellGID = graph(smPt, 0);
    const auto uIndex  = cellGID*1;

    const auto left0  = graph(smPt, 1);
    const auto front0 = graph(smPt, 2);
    const auto right0 = graph(smPt, 3);
    const auto back0  = graph(smPt, 4);

    if (left0 == -1){
      m_ghostLeft(gRow, 0) = -m_state(uIndex);
    }

    if (front0 == -1){
      m_ghostFront(gRow, 0) = -m_state(uIndex);
    }

    if (right0 == -1){
      m_ghostRight(gRow, 0) = -m_state(uIndex);
    }

    if (back0 == -1){
      m_ghostBack(gRow, 0) = -m_state(uIndex);
    }
  }

private:
  const int m_stencilSize;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  ghost_t & m_ghostLeft;
  ghost_t & m_ghostFront;
  ghost_t & m_ghostRight;
  ghost_t & m_ghostBack;
};

}}
#endif
