/*
//@HEADER
// ************************************************************************
//
// schemes_info.hpp
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

#ifndef PRESSIODEMOAPPS_RECONSTRUCTION_HPP_
#define PRESSIODEMOAPPS_RECONSTRUCTION_HPP_

namespace pressiodemoapps{

// here we have an enum also for a generic function
// reconstruction because this is not necessarily
// tied to inviscid/viscous stuff
enum class ReconstructionScheme{
  FirstOrder, Weno3, Weno5
};

enum class InviscidFluxReconstruction{
  FirstOrder, Weno3, Weno5
};
enum class InviscidFluxScheme{ Rusanov };
enum class ViscousFluxReconstruction{ FirstOrder };
enum class ViscousFluxScheme{ Central };


int reconstructionTypeToStencilSize(InviscidFluxReconstruction enIn){
  switch(enIn){
  case InviscidFluxReconstruction::FirstOrder: return 3;
  case InviscidFluxReconstruction::Weno3:	 return 5;
  case InviscidFluxReconstruction::Weno5:	 return 7;
  }
  return 0;
}

ReconstructionScheme toReconstructionScheme(InviscidFluxReconstruction enIn){
  switch(enIn){
  case InviscidFluxReconstruction::FirstOrder: return ReconstructionScheme::FirstOrder;
  case InviscidFluxReconstruction::Weno3:      return ReconstructionScheme::Weno3;
  case InviscidFluxReconstruction::Weno5:      return ReconstructionScheme::Weno5;
  }
  return {};
}

int reconstructionTypeToStencilSize(ViscousFluxReconstruction enIn){
  switch(enIn){
  case ViscousFluxReconstruction::FirstOrder: return 3;
  }
  return 0;
}

bool stencilSizeCompatibleWithInviscidFluxReconstruction(InviscidFluxReconstruction enIn,
							 int stencilSize){
  switch(enIn){
  case InviscidFluxReconstruction::FirstOrder: return stencilSize >= 3;
  case InviscidFluxReconstruction::Weno3:	 return stencilSize >= 5;
  case InviscidFluxReconstruction::Weno5:	 return stencilSize >= 7;
  }
  return false;
}

bool stencilSizeCompatibleWithViscousFluxReconstruction(ViscousFluxReconstruction enIn,
							int stencilSize){
  switch(enIn){
  case ViscousFluxReconstruction::FirstOrder: return stencilSize >= 3;
  }
  return false;
}

enum class FacePosition{
  Left, Front, Right, Back, Bottom, Top
};

enum class GradFdMode{
  ForwardTwoPt, BackwardTwoPt, CenterThreePt,
  ForwardThreePt, BackwardThreePt, CenterFivePt
};

enum class BoundaryFacesGradientScheme{
  OneSidedFdAutoStencil, // this selects the stencil cells automatically using the connectivity
  LSQAutoStencil
};

}//end namespace pressiodemoapps

#include "./weno.hpp"

#endif
