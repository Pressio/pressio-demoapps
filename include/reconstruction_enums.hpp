
#ifndef PRESSIODEMOAPPS_RECONSTRUCION_ENUMS_HPP_
#define PRESSIODEMOAPPS_RECONSTRUCION_ENUMS_HPP_

namespace pressiodemoapps{

enum class reconstructionEnum{
  firstOrder,
  fifthOrderWeno
};

int reconstructionEnumToStencilSize(reconstructionEnum enIn)
{
  switch(enIn)
    {
    case reconstructionEnum::firstOrder:     return 3;
    case reconstructionEnum::fifthOrderWeno: return 7;
    }
    return 0;
}

}
#endif
