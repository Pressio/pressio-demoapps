
#ifndef PRESSIODEMOAPPS_ADVECTION_1D_HPP_
#define PRESSIODEMOAPPS_ADVECTION_1D_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./adapter_mixins.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// enums identifying the problems
// ----------------------------------------------------------
enum class Advection1d{
  PeriodicLinear
};

}//end namespace pressiodemoapps

// this include is here because needs visiblity of the enums above
#include "./impl/advection_1d_prob_class.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// default problem
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<impladvection1d::EigenApp<mesh_t>>
  >
RetType
// bindings need unique naming or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
create_linear_advection_1d_problem_default_for_py
#else
create_problem_eigen
#endif
(const mesh_t & meshObj,
 Advection1d problemEnum,
 InviscidFluxReconstruction inviscidFluxRecEnum)
{

  using scalar_t = typename mesh_t::scalar_t;
  if (problemEnum == ::pressiodemoapps::Advection1d::PeriodicLinear){
    return RetType(impladvection1d::TagLinearAdvection{},
		   meshObj, inviscidFluxRecEnum,
		   InviscidFluxScheme::Rusanov,
		   static_cast<scalar_t>(1));
  }
  else{
    throw std::runtime_error("advection: invalid problem enum");
  }
}

// ----------------------------------------------------------
// problem with custom velocity
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<impladvection1d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
create_linear_advection_1d_problem_ov1_for_py
#else
create_linear_advection_1d_problem_eigen
#endif
(const mesh_t & meshObj,
 InviscidFluxReconstruction inviscidFluxRecEnum,
 InviscidFluxScheme inviscidFluxScheme,
 typename mesh_t::scalar_t velocity)
{

  return RetType(impladvection1d::TagLinearAdvection{},
		 meshObj, inviscidFluxRecEnum,
		 inviscidFluxScheme, velocity);
}

template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<impladvection1d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
create_linear_advection_1d_problem_ov2_for_py
#else
create_linear_advection_1d_problem_eigen
#endif
(const mesh_t & meshObj,
 InviscidFluxReconstruction inviscidFluxRecEnum,
 typename mesh_t::scalar_t velocity)
{
  return RetType(impladvection1d::TagLinearAdvection{},
		 meshObj, inviscidFluxRecEnum,
		 InviscidFluxScheme::Rusanov, velocity);
}

}//end namespace pressiodemoapps
#endif
