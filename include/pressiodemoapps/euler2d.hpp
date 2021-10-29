
#ifndef PRESSIODEMOAPPS_EULER2D_INC_HPP_
#define PRESSIODEMOAPPS_EULER2D_INC_HPP_

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
enum class Euler2d{
  PeriodicSmooth,
  KelvinHelmholtz,
  SedovFull,
  SedovSymmetry,
  Riemann,
  NormalShock,
  DoubleMachReflection,
  testingonlyneumann
};
}//end namespace pressiodemoapps

#include "./impl_apps/euler2d/euler2d_app_impl.hpp"

namespace pressiodemoapps{
namespace impl{
template<class mesh_t, class T>
T createEe2dImpl(const mesh_t & meshObj,
		 ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		 ::pressiodemoapps::Euler2d probEnum,
		 ::pressiodemoapps::InviscidFluxScheme fluxEnum,
		 const int icId)
{

  const auto stencilSize = meshObj.stencilSize();
  const auto check1 = stencilSizeCompatibleWithInviscidFluxReconstruction(recEnum, stencilSize);
  if (!check1){
    throw std::runtime_error
      ("Stencil size in the mesh object not compatible with desired inviscid flux reconstruction.");
  }

  if (probEnum == ::pressiodemoapps::Euler2d::DoubleMachReflection and
      recEnum == ::pressiodemoapps::InviscidFluxReconstruction::Weno5)
  {
    throw std::runtime_error
      ("Double mach reflection does NOT currently suppot Weno5.");
  }


  if (probEnum == ::pressiodemoapps::Euler2d::PeriodicSmooth)
  {
    if (!meshObj.isPeriodic()){
      throw std::runtime_error
      ("For euler2d::PeriodicSmooth, mesh must be periodic.");
    }
  }

  if (probEnum == ::pressiodemoapps::Euler2d::KelvinHelmholtz)
  {
    if (!meshObj.isPeriodic()){
      throw std::runtime_error
      ("For euler2d::KelvinHelmholtz, mesh must be periodic.");
    }
  }


  return T(meshObj, probEnum, recEnum, fluxEnum, icId);
}
} //end pressiodemoapps::impl

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler2d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum,
			::pressiodemoapps::InviscidFluxScheme fluxEnum,
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
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler2d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum,
			const int initCondIdentifier)
{
  return createProblemEigen(meshObj, probEnum, recEnum,
			    InviscidFluxScheme::Rusanov, initCondIdentifier);
}

template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler2d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return createProblemEigen(meshObj, probEnum, recEnum, InviscidFluxScheme::Rusanov, 1);
}
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createEuler2dForPyA(const mesh_t & meshObj,
		     ::pressiodemoapps::Euler2d probEnum,
		     ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		     ::pressiodemoapps::InviscidFluxScheme fluxEnum,
		     const int initCondIdentifier = 1)
{
  return impl::createEe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					 fluxEnum, initCondIdentifier);
}

template<class mesh_t, class T>
T createEuler2dForPyB(const mesh_t & meshObj,
		     ::pressiodemoapps::Euler2d probEnum,
		     ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		     const int icId = 1)
{
  return impl::createEe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					 InviscidFluxScheme::Rusanov, icId);
}


template<class mesh_t, class T>
T createEuler2dForPyC(const mesh_t & meshObj,
		     ::pressiodemoapps::Euler2d probEnum,
		     ::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return impl::createEe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					 InviscidFluxScheme::Rusanov, 1);
}
#endif

}
#endif
