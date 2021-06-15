
#ifndef PRESSIODEMOAPPS_RECONSTRUCTION_HPP_
#define PRESSIODEMOAPPS_RECONSTRUCTION_HPP_

namespace pressiodemoapps{

enum class ReconstructionType{
  firstOrder,
  fifthOrderWeno
};

int reconstructionTypeToStencilSize(ReconstructionType enIn)
{
  switch(enIn)
    {
    case ReconstructionType::firstOrder:     return 3;
    case ReconstructionType::fifthOrderWeno: return 7;
    }
    return 0;
}

#include "./weno.hpp"

}//end namespace pressiodemoapps
#endif
