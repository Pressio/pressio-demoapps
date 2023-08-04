/*
//@HEADER
// ************************************************************************
//
// euler_2d_prob_class.hpp
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

#ifndef PRESSIODEMOAPPS_FILL_GHOSTS_CUSTOM_BCS_HPP_
#define PRESSIODEMOAPPS_FILL_GHOSTS_CUSTOM_BCS_HPP_

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#include <omp.h>
#endif

namespace pressiodemoapps{
namespace impl{

template<class IndexType, class MeshType, class BCsFunctors, class JacFactors>
void fillJacFactorsCustomBCs(IndexType graphRow, int axis,
			     const MeshType & meshObj,
			     const BCsFunctors & bcFuncsHolder,
			     JacFactors & cellJacFactors,
			     int numDofPerCell)
{
  assert(axis <= 2);

  const auto & x = meshObj.get().viewX();
  const auto & y = meshObj.get().viewY();
  const auto & graph = meshObj.get().graph();
  auto currentCellGraphRow = graph.row(graphRow);
  const int cellGID = currentCellGraphRow[0];
  const auto myX = x(cellGID);
  const auto myY = y(cellGID);

  if (axis==1){
    if (meshObj.get().hasBdLeft2d(graphRow)){
      bcFuncsHolder(impl::GhostRelativeLocation::Left,
		      currentCellGraphRow, myX, myY, numDofPerCell,
		      cellJacFactors);
    }
    if (meshObj.get().hasBdRight2d(graphRow)){
      bcFuncsHolder(impl::GhostRelativeLocation::Right,
		      currentCellGraphRow, myX, myY, numDofPerCell,
		      cellJacFactors);
    }
  }
  else{
    if (meshObj.get().hasBdBack2d(graphRow)){
      bcFuncsHolder(impl::GhostRelativeLocation::Back,
		      currentCellGraphRow, myX, myY, numDofPerCell,
		      cellJacFactors);
    }

    if (meshObj.get().hasBdFront2d(graphRow)){
      bcFuncsHolder(impl::GhostRelativeLocation::Front,
		      currentCellGraphRow, myX, myY, numDofPerCell,
		      cellJacFactors);
    }
  }
}

template<
  class StateType, class ScalarType, class MeshType, class BCsFunctors,
  class GhostsLeft, class GhostsFront, class GhostsRight, class GhostsBack
  >
void fillGhostsUseCustomFunctors(const StateType & U,
				 const ScalarType currentTime,
				 const MeshType & meshObj,
				 const BCsFunctors & bcFuncsHolder,
				 GhostsLeft & ghostLeft,
				 GhostsFront & ghostFront,
				 GhostsRight & ghostRight,
				 GhostsBack & ghostBack,
				 int numDofPerCell)
{

  const auto & x = meshObj.get().viewX();
  const auto & y = meshObj.get().viewY();
  const auto & graph = meshObj.get().graph();
  const auto & rowsBd = meshObj.get().graphRowsOfCellsNearBd();

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
  for (decltype(rowsBd.size()) it=0; it<rowsBd.size(); ++it)
    {
      auto currentCellGraphRow = graph.row(rowsBd[it]);
      const int cellGID = rowsBd[it];
      const auto myX = x(cellGID);
      const auto myY = y(cellGID);

      /* IMPORTANT: keep the following as separate ifs wihtout ORs
	 because some cells might has ghosts on multiple sides so
	 we need these ifs not exclusive
      */
      if (meshObj.get().hasBdLeft2d(cellGID)){
	auto ghostVals = ghostLeft.row(it);
	bcFuncsHolder(impl::GhostRelativeLocation::Left,
		      it, currentCellGraphRow, myX, myY, U, numDofPerCell,
		      meshObj.get().dx(), ghostVals);
      }

      if (meshObj.get().hasBdFront2d(cellGID)){
	auto ghostVals = ghostFront.row(it);
	bcFuncsHolder(impl::GhostRelativeLocation::Front,
		      it, currentCellGraphRow, myX, myY, U, numDofPerCell,
		      meshObj.get().dy(), ghostVals);
      }

      if (meshObj.get().hasBdRight2d(cellGID)){
	auto ghostVals = ghostRight.row(it);
	bcFuncsHolder(impl::GhostRelativeLocation::Right,
		      it, currentCellGraphRow, myX, myY, U, numDofPerCell,
		      meshObj.get().dx(), ghostVals);
      }

      if (meshObj.get().hasBdBack2d(cellGID)){
	auto ghostVals = ghostBack.row(it);
	bcFuncsHolder(impl::GhostRelativeLocation::Back,
		      it, currentCellGraphRow, myX, myY, U, numDofPerCell,
		      meshObj.get().dy(), ghostVals);
      }
    }
}

}}
#endif
