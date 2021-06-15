
#ifndef PRESSIODEMOAPPS_EULER2D_INC_HPP_
#define PRESSIODEMOAPPS_EULER2D_INC_HPP_

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
#include "Eigen/Core"
#endif

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./reconstruction.hpp"
#include "flux_enum.hpp"

namespace pressiodemoapps{
enum class Euler2d{
  PeriodicSmooth,
  SedovFull,
  SedovSymmetry,
  Riemann,
  testingonlyneumann
};
}//end namespace pressiodemoapps

#include "./impl_apps/euler2d/euler2d_app_impl.hpp"

namespace pressiodemoapps{
namespace impl{
template<class mesh_t, class T>
T createEe2dImpl(const mesh_t & meshObj,
		 ::pressiodemoapps::ReconstructionType recEnum,
		 ::pressiodemoapps::Euler2d probEnum,
		 ::pressiodemoapps::FluxType fluxEnum,
		 const int icId)
{

  // // reconstruction order is specified, so we need to ensure
  // // the mesh object has stencil size that supports that
  // // e.g., firstOrder reconstruction can be done for 7-point stencil
  // // but 5-th order reconstruction cannot be done for 3-pt stencil
  // meshObj.checkStencilSupportsOrder(recEnum);

  if (probEnum == ::pressiodemoapps::Euler2d::PeriodicSmooth)
  {
    if (!meshObj.isPeriodic()){
      throw std::runtime_error
      ("For euler2d::PeriodicSmooth, mesh must be periodic.");
    }
  }

  return T(meshObj, probEnum, recEnum, fluxEnum, icId);
}
} //end pressiodemoapps::impl

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t>
auto create2dProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler2d probEnum,
			::pressiodemoapps::ReconstructionType recEnum,
			::pressiodemoapps::FluxType fluxEnum,
			int initCondIdentifier)
{

  using scalar_t = typename mesh_t::scalar_t;
  using p_t = ::pressiodemoapps::ee::impl::Euler2dAppT<
    scalar_t, mesh_t,
    Eigen::Matrix<scalar_t,-1,1>,
    Eigen::Matrix<scalar_t,-1,1>,
    Eigen::Matrix<scalar_t,-1,-1, Eigen::RowMajor>
    >;

  return impl::createEe2dImpl<mesh_t, p_t>(meshObj, recEnum, probEnum,
					   fluxEnum, initCondIdentifier);
}

template<class mesh_t>
auto create2dProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler2d probEnum,
			::pressiodemoapps::ReconstructionType recEnum,
			const int initCondIdentifier)
{
  return create2dProblemEigen(meshObj, probEnum, recEnum,
			    FluxType::Rusanov, initCondIdentifier);
}

template<class mesh_t>
auto create2dProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler2d probEnum,
			::pressiodemoapps::ReconstructionType recEnum)
{
  return create2dProblemEigen(meshObj, probEnum, recEnum, FluxType::Rusanov, 1);
}
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createEuler2dForPyA(const mesh_t & meshObj,
		     ::pressiodemoapps::Euler2d probEnum,
		     ::pressiodemoapps::ReconstructionType recEnum,
		     ::pressiodemoapps::FluxType fluxEnum,
		     const int initCondIdentifier = 1)
{
  return impl::createEe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					 fluxEnum, initCondIdentifier);
}

template<class mesh_t, class T>
T createEuler2dForPyB(const mesh_t & meshObj,
		     ::pressiodemoapps::Euler2d probEnum,
		     ::pressiodemoapps::ReconstructionType recEnum,
		     const int icId = 1)
{
  return impl::createEe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					 FluxType::Rusanov, icId);
}


template<class mesh_t, class T>
T createEuler2dForPyC(const mesh_t & meshObj,
		     ::pressiodemoapps::Euler2d probEnum,
		     ::pressiodemoapps::ReconstructionType recEnum)
{
  return impl::createEe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					 FluxType::Rusanov, 1);
}
#endif

}
#endif
