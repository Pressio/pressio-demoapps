
#ifndef PRESSIODEMOAPPS_EULER1D_INC_HPP_
#define PRESSIODEMOAPPS_EULER1D_INC_HPP_

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
#include "Eigen/Core"
#endif

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./reconstruction.hpp"
#include "./flux_enum.hpp"
#include "./euler_compute_energy.hpp"

namespace pressiodemoapps{
enum class Euler1d{
  PeriodicSmooth,
  Sod,
  Lax
};
}//end namespace pressiodemoapps

// this is here because it needs to see the enum above
#include "./impl_apps/euler1d/euler1d_app_impl.hpp"

namespace pressiodemoapps{

namespace impl{
template<class mesh_t, class T>
T createEuler1dImpl(const mesh_t & meshObj,
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

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createEuler1dForPyA(const mesh_t & meshObj,
		      ::pressiodemoapps::Euler1d probEnum,
		      ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		      ::pressiodemoapps::InviscidFluxScheme fluxEnum)
{
  return impl::createEuler1dImpl<mesh_t, T>(meshObj, probEnum, recEnum, fluxEnum);
}

template<class mesh_t, class T>
T createEuler1dForPyB(const mesh_t & meshObj,
		      ::pressiodemoapps::Euler1d probEnum,
		      ::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return impl::createEuler1dImpl<mesh_t, T>(meshObj, probEnum, recEnum,
					    InviscidFluxScheme::Rusanov);
}
#endif
} //end pressiodemoapps::impl


#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler1d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum,
			::pressiodemoapps::InviscidFluxScheme fluxEnum = InviscidFluxScheme::Rusanov)
{

  using scalar_t = typename mesh_t::scalar_t;
  //using p_t = ::pressiodemoapps::impl::Euler1dAppRhsOnly/*EigenEuler1dAppTWithJacobian*/<scalar_t, mesh_t>;
  using p_t = ::pressiodemoapps::impl::EigenEuler1dAppRhsOnly<scalar_t, mesh_t>;
  return impl::createEuler1dImpl<mesh_t, p_t>(meshObj, probEnum, recEnum, fluxEnum);
}
#endif

}//end namespace pressiodemoapps
#endif
