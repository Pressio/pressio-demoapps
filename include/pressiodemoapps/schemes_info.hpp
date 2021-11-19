
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

}//end namespace pressiodemoapps

#include "./weno.hpp"

#endif
