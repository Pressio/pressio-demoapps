
#ifndef PRESSIODEMOAPPS_EULER1D_INC_HPP_
#define PRESSIODEMOAPPS_EULER1D_INC_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./euler_compute_energy.hpp"
#include "./adapter_mixins.hpp"

namespace pressiodemoapps{
enum class Euler1d{
  PeriodicSmooth,
  Sod,
  Lax
};
}//end namespace pressiodemoapps

// this is here because it needs to see the enum above
#include "./impl/euler_1d_prob_class.hpp"

namespace pressiodemoapps{
namespace implee1d{

template<class mesh_t, class T>
T create1dImpl(const mesh_t & meshObj,
	       ::pressiodemoapps::Euler1d probEnum,
	       ::pressiodemoapps::InviscidFluxReconstruction recEnum,
	       ::pressiodemoapps::InviscidFluxScheme fluxEnum)
{

  const auto stencilSize = meshObj.stencilSize();
  const auto check1 = stencilSizeCompatibleWithInviscidFluxReconstruction(recEnum, stencilSize);
  if (!check1){
    throw std::runtime_error
      ("Stencil size in the mesh object not compatible with desired inviscid flux reconstruction.");
  }

  if (probEnum == ::pressiodemoapps::Euler1d::PeriodicSmooth){
    if (!meshObj.isPeriodic()){
      throw std::runtime_error("For periodicSmooth euler1d, mesh must be periodic.");
    }
  }

  return T(meshObj, probEnum, recEnum, fluxEnum);
}

#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createProblemForPyA(const mesh_t & meshObj,
		      ::pressiodemoapps::Euler1d probEnum,
		      ::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return implee1d::create1dImpl<mesh_t, T>(meshObj, probEnum,
					   recEnum, InviscidFluxScheme::Rusanov);
}

template<class mesh_t, class T>
T createProblemForPyB(const mesh_t & meshObj,
		      ::pressiodemoapps::Euler1d probEnum,
		      ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		      ::pressiodemoapps::InviscidFluxScheme fluxEnum)
{
  return implee1d::create1dImpl<mesh_t, T>(meshObj, probEnum, recEnum, fluxEnum);
}
#endif
} //end pressiodemoapps::impl

#if not defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler1d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum,
			::pressiodemoapps::InviscidFluxScheme fluxEnum = InviscidFluxScheme::Rusanov)
{
  using p_t = ::pressiodemoapps::implee1d::EigenEuler1dApp<mesh_t>;
  using RetType = PublicProblemMixinCpp<p_t>;
  return implee1d::create1dImpl<mesh_t, RetType>(meshObj, probEnum,
						 recEnum, fluxEnum);
}
#endif

}//end namespace pressiodemoapps
#endif
