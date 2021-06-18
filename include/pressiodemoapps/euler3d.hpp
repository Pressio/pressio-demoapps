
#ifndef PRESSIODEMOAPPS_EULER3D_INC_HPP_
#define PRESSIODEMOAPPS_EULER3D_INC_HPP_

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
#include "Eigen/Core"
#endif

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./reconstruction.hpp"
#include "flux_enum.hpp"

namespace pressiodemoapps{
enum class Euler3d{
  PeriodicSmooth,
  SedovSymmetry
};
}//end namespace pressiodemoapps

#include "./impl_apps/euler3d/euler3d_app_impl.hpp"

namespace pressiodemoapps{
namespace impl{
template<class mesh_t, class T>
T createEe3dImpl(const mesh_t & meshObj,
		 ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		 ::pressiodemoapps::Euler3d probEnum,
		 ::pressiodemoapps::InviscidFluxScheme fluxEnum,
		 const int icId)
{
  // // reconstruction order is specified, so we need to ensure
  // // the mesh object has stencil size that supports that
  // // e.g., FirstOrder reconstruction can be done for 7-point stencil
  // // but 5-th order reconstruction cannot be done for 3-pt stencil
  // meshObj.checkStencilSupportsOrder(recEnum);

  if (probEnum == ::pressiodemoapps::Euler3d::PeriodicSmooth)
  {
    if (!meshObj.isPeriodic()){
      throw std::runtime_error
      ("For euler3d::PeriodicSmooth, mesh must be periodic.");
    }
  }

  return T(meshObj, probEnum, recEnum, fluxEnum, icId);
}
} //end pressiodemoapps::impl

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler3d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum,
			::pressiodemoapps::InviscidFluxScheme fluxEnum,
			const int icId)
{
  using scalar_t = typename mesh_t::scalar_t;
  using p_t = ::pressiodemoapps::ee::impl::Euler3dAppT<
    scalar_t, mesh_t,
    Eigen::Matrix<scalar_t,-1,1>,
    Eigen::Matrix<scalar_t,-1,1>,
    Eigen::Matrix<scalar_t,-1,-1, Eigen::RowMajor>
    >;

  return impl::createEe3dImpl<mesh_t, p_t>(meshObj, recEnum, probEnum,
					   fluxEnum, icId);
}

template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler3d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return createProblemEigen(meshObj, probEnum, recEnum, InviscidFluxScheme::Rusanov, 1);
}
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createEuler3dForPyC(const mesh_t & meshObj,
		     ::pressiodemoapps::Euler3d probEnum,
		     ::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return impl::createEe3dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					 InviscidFluxScheme::Rusanov, 1);
}
#endif

}
#endif
