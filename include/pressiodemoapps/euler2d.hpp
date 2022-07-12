
#ifndef PRESSIODEMOAPPS_EULER2D_INC_HPP_
#define PRESSIODEMOAPPS_EULER2D_INC_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./euler_compute_energy.hpp"
#include "./adapter_cpp.hpp"
#include "./adapter_py.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// enums identifying the problems
// ----------------------------------------------------------
enum class Euler2d{
  PeriodicSmooth,
  KelvinHelmholtz,
  SedovFull,
  SedovSymmetry,
  Riemann,
  NormalShock,
  DoubleMachReflection,
  CrossShock,
  testingonlyneumann
};
}//end namespace pressiodemoapps

// this include is here because needs visiblity of the enums above
#include "./impl/euler_2d_prob_class.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// create a default problem
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impleuler2d::EigenApp<mesh_t>>
  >
RetType
// bindings need unique nameing or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  create_euler_2d_problem_default_for_py
#else
  create_problem_eigen
#endif
(const mesh_t & meshObj,
 Euler2d problemEnum,
 InviscidFluxReconstruction recEnum)
{
  return RetType(meshObj, problemEnum, recEnum,
		 InviscidFluxScheme::Rusanov, 1);
}

// ----------------------------------------------------------
// create a default problem with specific initial condition
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impleuler2d::EigenApp<mesh_t>>
  >
RetType
// bindings need unique nameing or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  create_euler_2d_problem_ov1_for_py
#else
  create_problem_eigen
#endif
(const mesh_t & meshObj,
 Euler2d problemEnum,
 InviscidFluxReconstruction recEnum,
 int icId)
{
  return RetType(meshObj, problemEnum, recEnum,
		 InviscidFluxScheme::Rusanov, icId);
}


// this crossshock one should really just be experimental
// and it is not documented on the website
template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impleuler2d::EigenApp<mesh_t>>
  >
RetType
// bindings need unique nameing or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  create_euler_2d_cross_shock_problem_for_py
#else
  create_cross_shock_problem_eigen
#endif
(const mesh_t & meshObj,
 InviscidFluxReconstruction recEnum,
 typename mesh_t::scalar_t density,
 typename mesh_t::scalar_t inletXVel,
 typename mesh_t::scalar_t bottomYVel)
{
  return RetType(impleuler2d::TagCrossShock{},
		 meshObj, recEnum,
		 InviscidFluxScheme::Rusanov, 1,
		 inletXVel, bottomYVel, density);
}

}//end namespace pressiodemoapps
#endif
