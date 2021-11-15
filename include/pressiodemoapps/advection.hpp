
#ifndef PRESSIODEMOAPPS_ADVECTION_HPP_
#define PRESSIODEMOAPPS_ADVECTION_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./adapter_mixins.hpp"

namespace pressiodemoapps{
enum class Advection1d{
  PeriodicLinear
};
}//end namespace pressiodemoapps

#include "./impl/advection_1d_prob_class.hpp"

namespace pressiodemoapps{
namespace impladv{

template<class mesh_t, class T>
T create1dImpl(const mesh_t & meshObj,
		     ::pressiodemoapps::Advection1d probEnum,
		     InviscidFluxReconstruction recEnum)
{
  const auto stencilSize = meshObj.stencilSize();
  const auto check1 = stencilSizeCompatibleWithInviscidFluxReconstruction(recEnum, stencilSize);
  if (!check1){
    throw std::runtime_error
      ("advection: stencil size in the mesh object not compatible with desired inviscid flux reconstruction.");
  }

  // the mesh should be periodic
  if (!meshObj.isPeriodic()){
    throw std::runtime_error("For periodic linear adv1d, mesh must be periodic.");
  }

  return T(meshObj, probEnum, recEnum);
}

#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<
  class mesh_t,
  class RetType = PublicProblemMixinPy<impladv::EigenAdvection1dApp<mesh_t>>
  >
RetType createProblemForPy(const mesh_t & meshObj,
			   ::pressiodemoapps::Advection1d probEnum,
			   ::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return impladv::create1dImpl<mesh_t, RetType>(meshObj, probEnum, recEnum);
}
#endif
} //end namespace impladv


#if not defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Advection1d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  using RetType = PublicProblemMixinCpp<impladv::EigenAdvection1dApp<mesh_t>>;
  return impladv::create1dImpl<mesh_t, RetType>(meshObj, probEnum, recEnum);
}
#endif

}//end namespace pressiodemoapps
#endif
