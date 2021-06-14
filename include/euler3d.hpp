
#ifndef PRESSIODEMOAPPS_EULER3D_INC_HPP_
#define PRESSIODEMOAPPS_EULER3D_INC_HPP_

#include <cmath>
#include <limits>

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
#include "Eigen/Core"
#endif

#include "./predicates/all.hpp"
#include "./resize.hpp"
#include "./extent.hpp"
#include "./weno.hpp"
#include "./reconstruction_enums.hpp"
#include "./mesh.hpp"
#include "./impl/eulerCommon/energy.hpp"
#include "./impl/eulerCommon/fluxes.hpp"
#include "./impl/stencil_filler.hpp"
#include "./impl/reconstructor_from_stencil.hpp"
#include "./impl/reconstructor_from_state.hpp"

namespace pressiodemoapps{
enum class euler3dproblemsEnum{
  PeriodicSmooth,
  SedovSymmetry
};

namespace impl{
template<class mesh_t, class T, typename ...Args>
T createEe3dImpl(const mesh_t & meshObj,
		 pressiodemoapps::reconstructionEnum enIn,
		 pressiodemoapps::euler3dproblemsEnum probEn,
		 Args && ... args)
{
  // reconstruction order is specified, so we need to ensure
  // the mesh object has stencil size that supports that
  // e.g., firstOrder reconstruction can be done for 7-point stencil
  // but 5-th order reconstruction cannot be done for 3-pt stencil
  meshObj.checkStencilSupportsOrder(enIn);

  if (probEn == pressiodemoapps::euler3dproblemsEnum::PeriodicSmooth)
  {
    if (!meshObj.isPeriodic()){
      throw std::runtime_error
      ("For euler3d::PeriodicSmooth, mesh must be periodic.");
    }
  }

  return T(meshObj, enIn, probEn, std::forward<Args>(args)...);
}
}} //end pressiodemoapps::impl


#include "./impl/euler3d/ghost_filler.hpp"
#include "./impl/euler3d/initial_condition.hpp"
#include "./impl/euler3d/euler3d_app_impl.hpp"

namespace pressiodemoapps{

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t, class scalar_t = double>
using Euler3dEigen =
  pressiodemoapps::ee::impl::Euler3dAppT<
  scalar_t, mesh_t,
  Eigen::Matrix<scalar_t,-1,1>,
  Eigen::Matrix<scalar_t,-1,1>,
  Eigen::Matrix<scalar_t,-1,-1, Eigen::RowMajor>
  >;

template<class mesh_t, typename ...Args>
auto createEuler3dEigen(const mesh_t & meshObj,
			pressiodemoapps::reconstructionEnum enIn,
			pressiodemoapps::euler3dproblemsEnum probEn,
			Args && ... args)
{
  using p_t = Euler3dEigen<mesh_t>;
  return impl::createEe3dImpl<mesh_t, p_t, Args...>(meshObj, enIn, probEn,
						    std::forward<Args>(args)...);

}
#endif

}
#endif
