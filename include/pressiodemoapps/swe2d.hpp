
#ifndef PRESSIODEMOAPPS_SWE2D_INC_HPP_
#define PRESSIODEMOAPPS_SWE2D_INC_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./adapter_mixins.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// enums identifying the problems
// ----------------------------------------------------------
enum class Swe2d{
  SlipWall,
};

}//end namespace pressiodemoapps

// this include is here because needs visiblity of the enums above
#include "./impl/swe_2d_prob_class.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// create default problem
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<implswe2d::EigenApp<mesh_t>>
  >
RetType
// bindings need unique naming or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  create_swe2d_problem_default_for_py
#else
  create_problem_eigen
#endif
(const mesh_t & meshObj,
 Swe2d problemEnum,
 InviscidFluxReconstruction inviscRecEn)
{

  using scalar_t = typename mesh_t::scalar_t;
  if (problemEnum == Swe2d::SlipWall)
  {
    return RetType(implswe2d::TagProblemSlipWall{},
		   meshObj, inviscRecEn,
		   ::pressiodemoapps::InviscidFluxScheme::Rusanov,
		   9.8,  // gravity
		   -3.0, // coriolis
		   0.125); // pulse mag
  }

  else{
    throw std::runtime_error("2D swe: invalid problem enum");
  }
}

// ----------------------------------------------------------
// custom coeffs
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<implswe2d::EigenApp<mesh_t>>
  >
RetType
// bindings need unique naming or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  create_slip_wall_swe2d_problem_ov1_for_py
#else
  create_slip_wall_swe_2d_problem_eigen
#endif
(const mesh_t & meshObj,
 InviscidFluxReconstruction inviscRecEn,
 typename mesh_t::scalar_t gravity,
 typename mesh_t::scalar_t coriolis,
 typename mesh_t::scalar_t pulseMagnitude)
{

  using scalar_t = typename mesh_t::scalar_t;
  return RetType(implswe2d::TagProblemSlipWall{},
		 meshObj, inviscRecEn,
		 ::pressiodemoapps::InviscidFluxScheme::Rusanov,
		 gravity, coriolis, pulseMagnitude);
}

}//end namespace pressiodemoapps
#endif
