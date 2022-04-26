
#ifndef PRESSIODEMOAPPS_ADVECTION_DIFFUSION_2D_HPP_
#define PRESSIODEMOAPPS_ADVECTION_DIFFUSION_2D_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./adapter_cpp.hpp"
#include "./adapter_py.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// enums identifying the problems
// ----------------------------------------------------------
enum class AdvectionDiffusion2d{
  Burgers
};

}//end namespace pressiodemoapps

// this include is here because needs visiblity of the enums above
#include "./impl/advection_diffusion_2d_prob_class.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// create default problem
// ----------------------------------------------------------
template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impladvdiff2d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// bindings need unique naming or we get error associated with overloads
create_advecdiffusion2d_problem_default_for_py
#else
create_problem_eigen
#endif
(const mesh_t & meshObj,
 AdvectionDiffusion2d problemEnum,
 InviscidFluxReconstruction inviscidFluxRecEnum,
 ViscousFluxReconstruction viscFluxRecEnum)
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

    return RetType(impladvdiff2d::TagBurgers{},
		   meshObj,
		   inviscidFluxRecEnum,
		   InviscidFluxScheme::Rusanov,
		   viscFluxRecEnum,
		   icPulseMagnitude, icSpread, diffusion,
		   icCenterX, icCenterY);
  }
  else{
    throw std::runtime_error("advection-diffusion2d: invalid problem enum");
  }
}

// ----------------------------------------------------------
// create problem with custom coefficients
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impladvdiff2d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// bindings need unique naming or we get error associated with overloads
create_burgers_2d_problem_ov1_for_py
#else
create_burgers_2d_problem_eigen
#endif
(const mesh_t & meshObj,
 InviscidFluxReconstruction inviscidFluxRecEnum,
 ViscousFluxReconstruction viscFluxRecEnum,
 typename mesh_t::scalar_t icPulseMagnitude,
 typename mesh_t::scalar_t icSpread,
 typename mesh_t::scalar_t diffusion,
 typename mesh_t::scalar_t icCenterX,
 typename mesh_t::scalar_t icCenterY)
{

  return RetType(impladvdiff2d::TagBurgers{},
		 meshObj,
		 inviscidFluxRecEnum,
		 InviscidFluxScheme::Rusanov,
		 viscFluxRecEnum,
		 icPulseMagnitude, icSpread, diffusion,
		 icCenterX, icCenterY);
}

} //end namespace pressiodemoapps
#endif
