
#ifndef PRESSIODEMOAPPS_ADVECTION_1D_HPP_
#define PRESSIODEMOAPPS_ADVECTION_1D_HPP_

/*
  1D advection enumerations and public APIs
*/

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./adapter_mixins.hpp"

namespace pressiodemoapps{

enum class Advection1d{
  PeriodicLinear
  // add more if needed
};

}//end namespace pressiodemoapps

// this include is here because needs visiblity of the enums above
#include "./impl/advection_1d_prob_class.hpp"

namespace pressiodemoapps{

#if not defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<impladv1d::EigenApp<mesh_t>>
  >
auto create_linear_advection1d_problem_eigen(const mesh_t & meshObj,
					     InviscidFluxReconstruction inviscidFluxRecEnum,
					     InviscidFluxScheme inviscidFluxScheme,
					     typename mesh_t::scalar_t velocity)
{

  return RetType(meshObj, inviscidFluxRecEnum, inviscidFluxScheme,
		 velocity, impladv1d::TagLinearAdvection{});
}

template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<impladv1d::EigenApp<mesh_t>>
  >
RetType create_problem_eigen(const mesh_t & meshObj,
			     Advection1d problemEnum,
			     InviscidFluxReconstruction inviscidFluxRecEnum,
			     InviscidFluxScheme inviscidFluxScheme = InviscidFluxScheme::Rusanov)
{

  using scalar_t = typename mesh_t::scalar_t;
  if (problemEnum == Advection1d::PeriodicLinear)
  {
    // default parameters
    return create_linear_advection1d_problem_eigen<mesh_t, RetType>(meshObj,
								    inviscidFluxRecEnum,
								    inviscidFluxScheme,
								    static_cast<scalar_t>(1));
  }
  else{
    throw std::runtime_error("advection: invalid problem enum");
  }
}
#endif

}//end namespace pressiodemoapps
#endif







// namespace impladv{

// template<class mesh_t, class T>
// T create1dImpl(const mesh_t & meshObj,
// 		     ::pressiodemoapps::Advection1d probEnum,
// 		     InviscidFluxReconstruction recEnum)
// {
//   const auto stencilSize = meshObj.stencilSize();
//   const auto check1 = stencilSizeCompatibleWithInviscidFluxReconstruction(recEnum, stencilSize);
//   if (!check1){
//     throw std::runtime_error
//       ("advection: stencil size in the mesh object not compatible with desired inviscid flux reconstruction.");
//   }

//   // the mesh should be periodic
//   if (!meshObj.isPeriodic()){
//     throw std::runtime_error("For periodic linear adv1d, mesh must be periodic.");
//   }

//   return T(meshObj, probEnum, recEnum);
// }

//#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// template<
//   class mesh_t,
//   class RetType = PublicProblemMixinPy<impladv::EigenAdvection1dApp<mesh_t>>
//   >
// RetType create_problem_for_py(const mesh_t & meshObj,
// 			   ::pressiodemoapps::Advection1d probEnum,
// 			   ::pressiodemoapps::InviscidFluxReconstruction recEnum)
// {
//   return impladv::create1dImpl<mesh_t, RetType>(meshObj, probEnum, recEnum);
// }
//#endif
//} //end namespace impladv
