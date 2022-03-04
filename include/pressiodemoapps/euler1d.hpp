
#ifndef PRESSIODEMOAPPS_EULER1D_INC_HPP_
#define PRESSIODEMOAPPS_EULER1D_INC_HPP_

/*
  Euler1d enumerations and public APIs
*/

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./euler_compute_energy.hpp"
#include "./adapter_mixins.hpp"

namespace pressiodemoapps{

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
// namespace implee1d{
// template<class mesh_t, class T>
// T create1dImpl(const mesh_t & meshObj,
// 	       ::pressiodemoapps::Euler1d probEnum,
// 	       ::pressiodemoapps::InviscidFluxReconstruction recEnum,
// 	       ::pressiodemoapps::InviscidFluxScheme fluxEnum)
// {

//   const auto stencilSize = meshObj.stencilSize();
//   const auto check1 = stencilSizeCompatibleWithInviscidFluxReconstruction(recEnum, stencilSize);
//   if (!check1){
//     throw std::runtime_error
//       ("Stencil size in the mesh object not compatible with desired inviscid flux reconstruction.");
//   }

//   if (probEnum == ::pressiodemoapps::Euler1d::PeriodicSmooth){
//     if (!meshObj.isPeriodic()){
//       throw std::runtime_error("For periodicSmooth euler1d, mesh must be periodic.");
//     }
//   }

//   return T(meshObj, probEnum, recEnum, fluxEnum);
// }

// #if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// template<class mesh_t, class T>
// T create_problem_for_pyA(const mesh_t & meshObj,
// 		      ::pressiodemoapps::Euler1d probEnum,
// 		      ::pressiodemoapps::InviscidFluxReconstruction recEnum)
// {
//   return implee1d::create1dImpl<mesh_t, T>(meshObj, probEnum,
// 					   recEnum, InviscidFluxScheme::Rusanov);
// }

// template<class mesh_t, class T>
// T create_problem_for_pyB(const mesh_t & meshObj,
// 		      ::pressiodemoapps::Euler1d probEnum,
// 		      ::pressiodemoapps::InviscidFluxReconstruction recEnum,
// 		      ::pressiodemoapps::InviscidFluxScheme fluxEnum)
// {
//   return implee1d::create1dImpl<mesh_t, T>(meshObj, probEnum, recEnum, fluxEnum);
// }
// #endif
//} //end pressiodemoapps::impl

#if not defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<impleuler1d::EigenApp<mesh_t>>
  >
RetType create_problem_eigen(const mesh_t & meshObj,
			     ::pressiodemoapps::Euler1d problemEnum,
			     ::pressiodemoapps::InviscidFluxReconstruction recEnum,
			     ::pressiodemoapps::InviscidFluxScheme fluxEnum = InviscidFluxScheme::Rusanov)
{
  return RetType(meshObj, problemEnum, recEnum, fluxEnum);
}
#endif

}//end namespace pressiodemoapps
#endif
