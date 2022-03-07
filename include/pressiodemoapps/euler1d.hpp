
#ifndef PRESSIODEMOAPPS_EULER1D_INC_HPP_
#define PRESSIODEMOAPPS_EULER1D_INC_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./euler_compute_energy.hpp"
#include "./adapter_mixins.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// enums identifying the problems
// ----------------------------------------------------------
enum class Euler1d{
  PeriodicSmooth,
  Sod,
  Lax,
  ShuOsher
};

}//end namespace pressiodemoapps

// this include is here because needs visiblity of the enums above
#include "./impl/euler_1d_prob_class.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// create a default problem
// ----------------------------------------------------------
template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<impleuler1d::EigenApp<mesh_t>>
  >
RetType
// bindings need unique nameing or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
create_euler1d_problem_ov1
#else
create_problem_eigen
#endif
(const mesh_t & meshObj,
 Euler1d problemEnum,
 InviscidFluxReconstruction recEnum)
{
  return RetType(meshObj, problemEnum, recEnum, InviscidFluxScheme::Rusanov);
}

}//end namespace pressiodemoapps
#endif
