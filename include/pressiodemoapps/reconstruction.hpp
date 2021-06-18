
#ifndef PRESSIODEMOAPPS_RECONSTRUCTION_HPP_
#define PRESSIODEMOAPPS_RECONSTRUCTION_HPP_

namespace pressiodemoapps{

enum class InviscidFluxReconstruction{
  FirstOrder,
  Weno5
};

int reconstructionTypeToStencilSize(InviscidFluxReconstruction enIn)
{
  switch(enIn)
    {
    case InviscidFluxReconstruction::FirstOrder:     return 3;
    case InviscidFluxReconstruction::Weno5: return 7;
    }
    return 0;
}

#include "./weno.hpp"

}//end namespace pressiodemoapps
#endif
