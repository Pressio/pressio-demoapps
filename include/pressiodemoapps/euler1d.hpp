
#ifndef PRESSIODEMOAPPS_EULER1D_INC_HPP_
#define PRESSIODEMOAPPS_EULER1D_INC_HPP_

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
#include "Eigen/Core"
#endif

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./reconstruction.hpp"
#include "flux_enum.hpp"

namespace pressiodemoapps{

enum class Euler1d{
  // for periodic, the user must provide a periodic domain and
  // a default smooth IC is provided
  PeriodicSmooth,

  // sod taken from: https://www.mdpi.com/2227-7390/6/10/211/pdf
  // x \in [-0.5, 0.5]
  // initial condition is provided by the app class
  Sod,

  // from https://www.sciencedirect.com/science/article/abs/pii/S0045793019300234
  // x \in [-5., 5.]
  // initial condition is provided by the app class
  Lax
};
}//end namespace pressiodemoapps

#include "./impl_apps/euler1d/euler1d_app_impl.hpp"

namespace pressiodemoapps{

namespace impl{
template<class mesh_t, class T>
T createEuler1dImpl(const mesh_t & meshObj,
		    ::pressiodemoapps::Euler1d probEnum,
		    ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		    ::pressiodemoapps::InviscidFluxScheme fluxEnum = InviscidFluxScheme::Rusanov)
{
  const auto stencilSize = meshObj.stencilSize();
  const auto check1 = stencilSizeCompatibleWithInviscidFluxReconstruction(recEnum, stencilSize);
  if (!check1){
    throw std::runtime_error
      ("Stencil size in the mesh object not compatible with desired inviscid flux reconstruction.");
  }

  if (probEnum == ::pressiodemoapps::Euler1d::PeriodicSmooth){
    if (!meshObj.isPeriodic()){
      throw std::runtime_error
      ("For periodicSmooth euler1d, mesh must be periodic.");
    }
  }

  return T(meshObj, probEnum, recEnum, fluxEnum);
}
} //end pressiodemoapps::impl

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler1d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum,
			::pressiodemoapps::InviscidFluxScheme fluxEnum = InviscidFluxScheme::Rusanov)
{

  using scalar_t = typename mesh_t::scalar_t;
  using p_t = ::pressiodemoapps::ee::impl::Euler1dAppT<
    scalar_t, mesh_t,
    Eigen::Matrix<scalar_t,-1,1>,
    Eigen::Matrix<scalar_t,-1,1>,
    Eigen::Matrix<scalar_t,-1,1>
    >;
  return impl::createEuler1dImpl<mesh_t, p_t>(meshObj, probEnum, recEnum, fluxEnum);
}
#endif

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
  return impl::createEuler1dImpl<mesh_t, T>(meshObj, probEnum,
					    recEnum, InviscidFluxScheme::Rusanov);
}
#endif

}//end namespace pressiodemoapps
#endif
