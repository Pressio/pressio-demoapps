
#ifndef PRESSIODEMOAPPS_ADVECTION_HPP_
#define PRESSIODEMOAPPS_ADVECTION_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./reconstruction.hpp"

namespace pressiodemoapps{
enum class Advection1d{
  PeriodicLinear
};
}//end namespace pressiodemoapps

#include "./impl_apps/advection1d/linear_adv_impl.hpp"

namespace pressiodemoapps{
namespace impl{
template<class mesh_t, class T>
T createLinAdv1dImpl(const mesh_t & meshObj,
		     ::pressiodemoapps::Advection1d probEnum,
		     InviscidFluxReconstruction recEnum)
{
  // // reconstruction order is specified, so we need to ensure
  // // the mesh object has stencil size that supports that
  // // e.g., FirstOrder reconstruction can be done for 7-point stencil
  // // but 5-th order reconstruction cannot be done for 3-pt stencil
  // meshObj.checkStencilSupportsOrder(recEnum);

  // the mesh should be periodic
  if (!meshObj.isPeriodic()){
    throw std::runtime_error
      ("For periodic linear adv1d, mesh must be periodic.");
  }

  return T(meshObj, probEnum, recEnum);
}
} //end impl

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Advection1d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  using scalar_t = typename mesh_t::scalar_t;
  using return_type = ::pressiodemoapps::ad::impl::AdvectionAppT<
    scalar_t, mesh_t, Eigen::Matrix<scalar_t,-1,1>, Eigen::Matrix<scalar_t,-1,1>
    >;
  return impl::createLinAdv1dImpl<mesh_t, return_type>(meshObj, probEnum, recEnum);
}
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createAdv1dForPy(const mesh_t & meshObj,
		   ::pressiodemoapps::Advection1d probEnum,
		   ::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return impl::createLinAdv1dImpl<mesh_t, T>(meshObj, probEnum, recEnum);
}
#endif

}//end namespace pressiodemoapps
#endif
