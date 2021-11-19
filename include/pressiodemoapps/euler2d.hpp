
#ifndef PRESSIODEMOAPPS_EULER2D_INC_HPP_
#define PRESSIODEMOAPPS_EULER2D_INC_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./euler_compute_energy.hpp"
#include "./adapter_mixins.hpp"

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

#include "./impl/euler_2d_prob_class.hpp"

namespace pressiodemoapps{
namespace implee2d{

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


#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createProblemForPyA(const mesh_t & meshObj,
		      ::pressiodemoapps::Euler2d probEnum,
		      ::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return implee2d::createEe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					     InviscidFluxScheme::Rusanov, 1);
}

template<class mesh_t, class T>
T createProblemForPyB(const mesh_t & meshObj,
		      ::pressiodemoapps::Euler2d probEnum,
		      ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		      const int ic)
{
  return implee2d::createEe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					     InviscidFluxScheme::Rusanov, ic);
}
#endif

} //end pressiodemoapps::impl

#if not defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler2d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum,
			::pressiodemoapps::InviscidFluxScheme fluxEnum,
			int initCondIdentifier)
{
  using p_t = ::pressiodemoapps::ee::impl::EigenEuler2dApp<mesh_t>;
  using RetType = PublicProblemMixinCpp<p_t>;
  return implee2d::createEe2dImpl<mesh_t, RetType>(meshObj, recEnum, probEnum,
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
  return createProblemEigen(meshObj, probEnum, recEnum,
			    InviscidFluxScheme::Rusanov, 1);
}
#endif

}//end namespace pressiodemoapps
#endif
