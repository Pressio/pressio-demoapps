
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
  RayleighTaylor,
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
T create_problem_for_pyA(const mesh_t & meshObj,
			 ::pressiodemoapps::Euler2d probEnum,
			 ::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return implee2d::createEe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					     InviscidFluxScheme::Rusanov, 1);
}

template<class mesh_t, class T>
T create_problem_for_pyB(const mesh_t & meshObj,
			 ::pressiodemoapps::Euler2d probEnum,
			 ::pressiodemoapps::InviscidFluxReconstruction recEnum,
			 const int ic)
{
  return implee2d::createEe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					     InviscidFluxScheme::Rusanov, ic);
}

template<class mesh_t, class T>
T create_problem_for_pyC(const mesh_t & meshObj,
			 ::pressiodemoapps::Euler2d probEnum,
			 ::pressiodemoapps::InviscidFluxReconstruction recEnum,
			 const typename mesh_t::scalar_t amplitude)
{
  return T(meshObj, probEnum, recEnum, 5./3., amplitude);
}
#endif

} //end pressiodemoapps::implee2d

#if not defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t>
auto create_problem_eigen(const mesh_t & meshObj,
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
auto create_problem_eigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler2d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum,
			const int initCondIdentifier)
{
  return create_problem_eigen(meshObj, probEnum, recEnum,
			    InviscidFluxScheme::Rusanov, initCondIdentifier);
}

template<class mesh_t>
auto create_problem_eigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler2d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return create_problem_eigen(meshObj, probEnum, recEnum,
			    InviscidFluxScheme::Rusanov, 1);
}

// for RTI problem
template<class mesh_t>
auto create_problem_eigen(const mesh_t & meshObj,
			  ::pressiodemoapps::Euler2d probEnum,
			  ::pressiodemoapps::InviscidFluxReconstruction recEnum,
			  typename mesh_t::scalar_t amplitude)
{
  if (probEnum != ::pressiodemoapps::Euler2d::RayleighTaylor){
    throw std::runtime_error("constructor valid only for RayleighTaylor");
  }

  using sc_t = typename mesh_t::scalar_t;
  const auto gamma = static_cast<sc_t>(5)/static_cast<sc_t>(3);
  using p_t = ::pressiodemoapps::ee::impl::EigenEuler2dApp<mesh_t>;
  using RetType = PublicProblemMixinCpp<p_t>;
  return RetType(meshObj, probEnum, recEnum,
		 gamma, amplitude);
}

template<class mesh_t>
auto create_problem_eigen(const mesh_t & meshObj,
			  ::pressiodemoapps::Euler2d probEnum,
			  ::pressiodemoapps::InviscidFluxReconstruction recEnum,
			  typename mesh_t::scalar_t gamma,
			  typename mesh_t::scalar_t amplitude)
{
  if (probEnum != ::pressiodemoapps::Euler2d::RayleighTaylor){
    throw std::runtime_error("constructor valid only for RayleighTaylor");
  }

  using sc_t = typename mesh_t::scalar_t;
  using p_t = ::pressiodemoapps::ee::impl::EigenEuler2dApp<mesh_t>;
  using RetType = PublicProblemMixinCpp<p_t>;
  return RetType(meshObj, probEnum, recEnum,
		 gamma, amplitude);
}
#endif

}//end namespace pressiodemoapps
#endif
