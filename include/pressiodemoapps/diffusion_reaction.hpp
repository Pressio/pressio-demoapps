
#ifndef PRESSIODEMOAPPS_DIFFUSION_REACTION_HPP_
#define PRESSIODEMOAPPS_DIFFUSION_REACTION_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./adapter_mixins.hpp"

namespace pressiodemoapps{
enum class DiffusionReaction1d{
  /*
    ds/dt = D d^2s/d^2x + k*s^2 + u(x, t)

    BC: ghost cells are set such that s is zero at boundary

    D, k, u(x, t) can be provided to constructor
    u(x, t) must be a functor:
      void operator()(const scalar_t & x,
		      const scalar_t & time,
		      scalar_t & value);

    Default:
      D,k = 0.01, 0.01
      u(x, t) = sin(M_PI*x) * x*x * 4.*std::cos(4.*M_PI*x);
   */
  ProblemA
};

// enum class DiffusionReaction2d{
//   /*
//     ds/dt = D (d^2s/d^2x + d^2s/d^2y) + k*s^2 + u(x, y, t)

//     D, k, u(x, y, t) can be provided to constructor
//     u(x, y, t) must be a functor:
//       void operator()(const scalar_t & x,
// 		      const scalar_t & y,
// 		      const scalar_t & time,
// 		      scalar_t & value);

//     BC: ghost cells are set such that s is zero at boundary

//     Default:
//     D, k = 0.01, 0.01
//     u(x, y, t) = std::sin(M_PI*x*(y-0.2)) * 4.*std::sin(4.*M_PI*y*x);

//    */
//   ProblemA
// };
}//end namespace pressiodemoapps

#include "./impl/diffusion_reaction_1d_prob_class.hpp"
//#include "./impl/diffusion_reaction_2d_prob_class.hpp"

namespace pressiodemoapps{
namespace impldiffreac{

template<class mesh_t>
void checkStencilAdmissibility(const mesh_t & meshObj,
			       ::pressiodemoapps::ViscousFluxReconstruction recEnum)
{
  const auto stencilSize = meshObj.stencilSize();
  const auto check1 = stencilSizeCompatibleWithViscousFluxReconstruction(recEnum, stencilSize);
  if (!check1){
    throw std::runtime_error
      ("Stencil size in the mesh object not compatible with desired viscous flux reconstruction.");
  }
}

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createProblemForPyA(const mesh_t & meshObj,
		      ::pressiodemoapps::DiffusionReaction1d probEnum,
		      ::pressiodemoapps::ViscousFluxReconstruction recEnum)
{
  impldiffreac::checkStencilAdmissibility(meshObj, recEnum);
  return T(meshObj, probEnum, recEnum);
}

template<class mesh_t, class T>
T createProblemForPyB(const mesh_t & meshObj,
		      ::pressiodemoapps::DiffusionReaction1d probEnum,
		      ::pressiodemoapps::ViscousFluxReconstruction recEnum,
		      typename mesh_t::scalar_t diffusionCoeff,
		      typename mesh_t::scalar_t reactionCoeff)
{
  impldiffreac::checkStencilAdmissibility(meshObj, recEnum);
  return T(meshObj, probEnum, recEnum, diffusionCoeff, reactionCoeff);
}

template<class mesh_t, class T>
T createProblemForPyC(const mesh_t & meshObj,
		      ::pressiodemoapps::DiffusionReaction1d probEnum,
		      ::pressiodemoapps::ViscousFluxReconstruction recEnum,
		      pybind11::object pyFunctor,
		      typename mesh_t::scalar_t diffusionCoeff,
		      typename mesh_t::scalar_t reactionCoeff)
{
  impldiffreac::checkStencilAdmissibility(meshObj, recEnum);
  using scalar_t = typename mesh_t::scalar_t;
  auto sourceWrapper = [=](const scalar_t & x, const scalar_t & evaltime, scalar_t & value){
    value = pyFunctor.attr("__call__")(x, evaltime).template cast<scalar_t>();
  };

  return T(meshObj, probEnum, recEnum, sourceWrapper, diffusionCoeff, reactionCoeff);
}
#endif
} //end namespace impldiffreac

#if not defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class ...Args>
auto createExplicitProblemEigen(const mesh_t & meshObj,
				::pressiodemoapps::DiffusionReaction1d probEnum,
				::pressiodemoapps::ViscousFluxReconstruction recEnum,
				Args && ... args)
{
  impldiffreac::checkStencilAdmissibility(meshObj, recEnum);
  using RetType = ExplicitAdapterMixin<impldiffreac::EigenDiffReac1dApp<mesh_t>>;
  return RetType(meshObj, probEnum, recEnum, std::forward<Args>(args)...);
}

template<class mesh_t, class ...Args>
auto createImplicitProblemEigen(const mesh_t & meshObj,
				::pressiodemoapps::DiffusionReaction1d probEnum,
				::pressiodemoapps::ViscousFluxReconstruction recEnum,
				Args && ... args)
{
  impldiffreac::checkStencilAdmissibility(meshObj, recEnum);
  using RetType = ImplicitAdapterMixinCpp<impldiffreac::EigenDiffReac1dApp<mesh_t>>;
  return RetType(meshObj, probEnum, recEnum, std::forward<Args>(args)...);
}
#endif

} //end namespace pressiodemoapps
#endif
