
#ifndef PRESSIODEMOAPPS_ADVECTION_DIFFUSION_2D_HPP_
#define PRESSIODEMOAPPS_ADVECTION_DIFFUSION_2D_HPP_

/*
  2D advection-diffusion enumerations and public APIs
*/

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./adapter_mixins.hpp"

namespace pressiodemoapps{

enum class AdvectionDiffusion2d{
  Burgers
  // add more if needed
};

}//end namespace pressiodemoapps

// this include is here because needs visiblity of the enums above
#include "./impl/advection_diffusion_2d_prob_class.hpp"

namespace pressiodemoapps{
namespace impladvdiff2d{

// #ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
// template<class mesh_t, class T, class ProbType>
// T create_problem_for_pyA(const mesh_t & meshObj,
// 			 ProbType probEnum,
// 			 ::pressiodemoapps::InviscidFluxReconstruction invFluxRecEnum,
// 			 typename mesh_t::scalar_t icPulseMagnitude,
// 			 typename mesh_t::scalar_t icPulseSpread,
// 			 typename mesh_t::scalar_t diffusionCoeff,
// 			 typename mesh_t::scalar_t x0,
// 			 typename mesh_t::scalar_t y0)
// {
//   return T(meshObj, probEnum, invFluxRecEnum,
// 	   ::pressiodemoapps::ViscousFluxReconstruction::FirstOrder,
// 	   icPulseMagnitude, icPulseSpread, diffusionCoeff, x0, y0);
// }
// #endif

} //end namespace impladvdiff2d


//
// public API for constructing problems
//
#if not defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<impladvdiff2d::EigenApp<mesh_t>>
  >
auto create_burgers2d_problem_eigen(const mesh_t & meshObj,
				    ::pressiodemoapps::InviscidFluxReconstruction inviscidFluxRecEnum,
				    ::pressiodemoapps::InviscidFluxScheme inviscidFluxScheme,
				    ::pressiodemoapps::ViscousFluxReconstruction viscFluxRecEnum,
				    typename mesh_t::scalar_t icPulseMagnitude,
				    typename mesh_t::scalar_t icSpread,
				    typename mesh_t::scalar_t diffusion,
				    typename mesh_t::scalar_t icCenterX,
				    typename mesh_t::scalar_t icCenterY)
{

  return RetType(meshObj,
		 inviscidFluxRecEnum,
		 inviscidFluxScheme,
		 viscFluxRecEnum,
		 icPulseMagnitude, icSpread, diffusion,
		 icCenterX, icCenterY,
		 impladvdiff2d::TagBurgers{});
}

template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<impladvdiff2d::EigenApp<mesh_t>>
  >
RetType create_problem_eigen(const mesh_t & meshObj,
			     ::pressiodemoapps::AdvectionDiffusion2d problemEnum,
			     ::pressiodemoapps::InviscidFluxReconstruction inviscidFluxRecEnum,
			     ::pressiodemoapps::ViscousFluxReconstruction viscFluxRecEnum)
{

  using scalar_t = typename mesh_t::scalar_t;
  if (problemEnum == ::pressiodemoapps::AdvectionDiffusion2d::Burgers)
  {
    // default parameters for burgers2d:
    const auto icPulseMagnitude = static_cast<scalar_t>(0.5);
    const auto icSpread = static_cast<scalar_t>(0.15);
    const auto diffusion = static_cast<scalar_t>(0.00001);
    const auto icCenterX = static_cast<scalar_t>(0.0);
    const auto icCenterY = static_cast<scalar_t>(-0.2);

    return create_burgers2d_problem_eigen<mesh_t, RetType>(meshObj,
							   inviscidFluxRecEnum,
							   ::pressiodemoapps::InviscidFluxScheme::Rusanov,
							   viscFluxRecEnum,
							   icPulseMagnitude,
							   icSpread, diffusion,
							   icCenterX, icCenterY);
  }
  else{
    throw std::runtime_error("advection-diffusion: invalid problem enum");
  }
}
#endif

} //end namespace pressiodemoapps
#endif
