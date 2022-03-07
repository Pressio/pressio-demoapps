
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
// create problem with custom velocity
// ----------------------------------------------------------
template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<impladvection1d::EigenApp<mesh_t>>
  >
auto create_linear_advection_1d_problem_eigen(const mesh_t & meshObj,
					      InviscidFluxReconstruction inviscidFluxRecEnum,
					      InviscidFluxScheme inviscidFluxScheme,
					      typename mesh_t::scalar_t velocity)
{

  return RetType(impladvection1d::TagLinearAdvection{},
		 meshObj, inviscidFluxRecEnum, inviscidFluxScheme, velocity);
}

// ----------------------------------------------------------
// create default problem
// ----------------------------------------------------------
template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<impladvection1d::EigenApp<mesh_t>>
  >
RetType
// bindings need unique naming or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
create_linadv1d_problem_ov1
#else
create_problem_eigen
#endif
(const mesh_t & meshObj,
 Advection1d problemEnum,
 InviscidFluxReconstruction inviscidFluxRecEnum)
{

  using scalar_t = typename mesh_t::scalar_t;
  if (problemEnum == ::pressiodemoapps::Advection1d::PeriodicLinear)
    {
      // default velocity is 1.
      return create_linear_advection_1d_problem_eigen<
	mesh_t, RetType>(meshObj, inviscidFluxRecEnum,
			 ::pressiodemoapps::InviscidFluxScheme::Rusanov,
			 static_cast<scalar_t>(1));
    }
  else{
    throw std::runtime_error("advection: invalid problem enum");
  }
}

}//end namespace pressiodemoapps
#endif
