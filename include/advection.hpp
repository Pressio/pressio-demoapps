
#ifndef PRESSIODEMOAPPS_ADVECTION_HPP_
#define PRESSIODEMOAPPS_ADVECTION_HPP_

#include "./predicates/all.hpp"
#include "./resize.hpp"
#include "./extent.hpp"
#include "./weno.hpp"
#include "./reconstruction_enums.hpp"
#include "./mesh.hpp"
#include "./impl/stencil_filler.hpp"
#include "./impl/reconstructor_from_stencil.hpp"
#include "./impl/reconstructor_from_state.hpp"
#include "./impl/advection1d/fluxes.hpp"
#include "./impl/advection1d/linear_adv_impl.hpp"

namespace pressiodemoapps{

namespace impl{
template<class mesh_t, class T>
T createLinAdv1dImpl(const mesh_t & meshObj,
		     pressiodemoapps::reconstructionEnum enIn)
{
  // reconstruction order is specified, so we need to ensure
  // the mesh object has stencil size that supports that
  // e.g., firstOrder reconstruction can be done for 7-point stencil
  // but 5-th order reconstruction cannot be done for 3-pt stencil
  meshObj.checkStencilSupportsOrder(enIn);

  // the mesh should be periodic
  if (!meshObj.isPeriodic()){
    throw std::runtime_error
      ("For periodic linear adv1d, mesh must be periodic.");
  }

  return T(meshObj, enIn);
}
} //end impl

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t, class scalar_t = double>
using PeriodicLinearAdvection1dEigen =
  pressiodemoapps::ad::impl::LinearAdvT<
  scalar_t, mesh_t,
  Eigen::Matrix<scalar_t,-1,1>,
  Eigen::Matrix<scalar_t,-1,1>
  >;

template<class mesh_t>
auto createPeriodicLinearAdvection1dEigen(const mesh_t & meshObj,
					  pressiodemoapps::reconstructionEnum enIn)
{
  using p_t = PeriodicLinearAdvection1dEigen<mesh_t>;
  return impl::createLinAdv1dImpl<mesh_t, p_t>(meshObj, enIn);
}
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createPeriodicLinearAdvection1dForPy(const mesh_t & meshObj,
					  pressiodemoapps::reconstructionEnum enIn)
{
  return impl::createLinAdv1dImpl<mesh_t, T>(meshObj, enIn);
}
#endif

}//end namespace pressiodemoapps
#endif
