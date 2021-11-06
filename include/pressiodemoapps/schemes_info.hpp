
#ifndef PRESSIODEMOAPPS_RECONSTRUCTION_HPP_
#define PRESSIODEMOAPPS_RECONSTRUCTION_HPP_

namespace pressiodemoapps{

enum class InviscidFluxScheme{
  Rusanov
};

enum class InviscidFluxReconstruction{
  FirstOrder,
  Weno3,
  Weno5
};

enum class ViscousFluxScheme{
  Central
};

enum class ViscousFluxReconstruction{
  FirstOrder
};

int reconstructionTypeToStencilSize(InviscidFluxReconstruction enIn)
{
  switch(enIn)
    {
    case InviscidFluxReconstruction::FirstOrder: return 3;
    case InviscidFluxReconstruction::Weno3:	 return 5;
    case InviscidFluxReconstruction::Weno5:	 return 7;
    }
    return 0;
}

int reconstructionTypeToStencilSize(ViscousFluxReconstruction enIn)
{
  switch(enIn)
    {
    case ViscousFluxReconstruction::FirstOrder: return 3;
    }
    return 0;
}

bool stencilSizeCompatibleWithInviscidFluxReconstruction
(InviscidFluxReconstruction enIn, int stencilSize)
{
  switch(enIn)
    {
    case InviscidFluxReconstruction::FirstOrder: return stencilSize >= 3;
    case InviscidFluxReconstruction::Weno3:	 return stencilSize >= 5;
    case InviscidFluxReconstruction::Weno5:	 return stencilSize >= 7;
    }
    return false;
}

bool stencilSizeCompatibleWithViscousFluxReconstruction
(ViscousFluxReconstruction enIn, int stencilSize)
{
  switch(enIn)
    {
    case ViscousFluxReconstruction::FirstOrder: return stencilSize >= 3;
    }
  return false;
}

#include "./weno.hpp"

}//end namespace pressiodemoapps
#endif
