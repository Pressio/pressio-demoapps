
#ifndef PRESSIODEMOAPPS_EULER2D_INC_HPP_
#define PRESSIODEMOAPPS_EULER2D_INC_HPP_

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
enum class euler2dproblemsEnum{
  PeriodicSmooth,
  SedovFull,
  SedovSymmetry,
  Riemann,
  testingonlyneumann
};

namespace impl{
template<class mesh_t, class T, typename ...Args>
T createEe2dImpl(const mesh_t & meshObj,
		 pressiodemoapps::reconstructionEnum enIn,
		 pressiodemoapps::euler2dproblemsEnum probEn,
		 Args && ... args)
{
  // reconstruction order is specified, so we need to ensure
  // the mesh object has stencil size that supports that
  // e.g., firstOrder reconstruction can be done for 7-point stencil
  // but 5-th order reconstruction cannot be done for 3-pt stencil
  meshObj.checkStencilSupportsOrder(enIn);

  if (probEn == pressiodemoapps::euler2dproblemsEnum::PeriodicSmooth)
  {
    if (!meshObj.isPeriodic()){
      throw std::runtime_error
      ("For euler2d::PeriodicSmooth, mesh must be periodic.");
    }
  }

  return T(meshObj, enIn, probEn, std::forward<Args>(args)...);
}
}} //end pressiodemoapps::impl


#include "./impl/euler2d/euler2d_app_impl.hpp"

namespace pressiodemoapps{

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t, class scalar_t = double>
using Euler2dEigen =
  pressiodemoapps::ee::impl::Euler2dAppT<
  scalar_t, mesh_t,
  Eigen::Matrix<scalar_t,-1,1>,
  Eigen::Matrix<scalar_t,-1,1>,
  Eigen::Matrix<scalar_t,-1,-1, Eigen::RowMajor>
  >;

template<class mesh_t, typename ...Args>
auto createEuler2dEigen(const mesh_t & meshObj,
			pressiodemoapps::reconstructionEnum enIn,
			pressiodemoapps::euler2dproblemsEnum probEn,
			Args && ... args)
{
  using p_t = Euler2dEigen<mesh_t>;
  return impl::createEe2dImpl<mesh_t, p_t, Args...>(meshObj, enIn, probEn,
						    std::forward<Args>(args)...);

}
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createEuler2dForPy(const mesh_t & meshObj,
		     pressiodemoapps::reconstructionEnum enIn,
		     pressiodemoapps::euler2dproblemsEnum probEn,
		     const int initCondIdentifier)
{
  return impl::createEe2dImpl<mesh_t, T>(meshObj, enIn,
					 probEn, initCondIdentifier);
}
#endif

}
#endif
