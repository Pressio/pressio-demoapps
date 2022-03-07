
#ifndef PRESSIODEMOAPPS_DIFFUSION_REACTION_HPP_
#define PRESSIODEMOAPPS_DIFFUSION_REACTION_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./adapter_mixins.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// enums identifying the problems
// ----------------------------------------------------------
enum class DiffusionReaction2d
{
  ProblemA,
  /*
    ds/dt = D (d^2s/d^2x + d^2s/d^2y) + k*s^2 + u(x, y, t)

    BC: ghost cells are set such that s is zero at boundary

    D, k, u(x, y, t) can be provided to constructor
    u(x, y, t) must be a functor:
      void operator()(const scalar_t & x,
		      const scalar_t & y,
		      const scalar_t & time,
		      scalar_t & value);

    Default:
    D, k = 0.01, 0.01
    u(x, y, t) = std::sin(M_PI*x*(y-0.2)) * 4.*std::sin(4.*M_PI*y*x);

   */

  GrayScott
  /*
    two species:
    da/dt = Da Laplacian(a) + F(1-a) - a*b^2
    db/dt = Db Laplacian(v) + (F+k)b + a*b^2

    F = feeding rate, k = kill rate
    Da, Db = diffusion coeffs
    periodic BC
   */
};

// ----------------------------------------------------------
// default source terms
// ----------------------------------------------------------
namespace impldiffusionreaction2d{
template<class scalar_type>
struct DefaultSourceProblemA
{
  void operator()(const scalar_type & x,
		  const scalar_type & y,
		  const scalar_type & evaltime,
		  scalar_type & value)
  {
    (void) evaltime;
    value = std::sin(M_PI*x*(y-0.2)) * 4.*std::sin(4.*M_PI*y*x);
  }
};
}//end namespace impldiffusionreaction1d
}//end namespace pressiodemoapps


// this include is here because needs visiblity of the enums above
#include "./impl/diffusion_reaction_2d_prob_class.hpp"


namespace pressiodemoapps{

// ----------------------------------------------------------
// create default problem
// ----------------------------------------------------------
template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<impldiffusionreaction2d::EigenApp<mesh_t>>
  >
RetType
// bindings need unique naming or we get error associated with overloads
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
create_diffreac2d_problem_ov1
#else
create_problem_eigen
#endif
(const mesh_t & meshObj,
 DiffusionReaction2d problemEnum,
 ViscousFluxReconstruction viscRecEn)
{

  using scalar_t = typename mesh_t::scalar_t;
  if (problemEnum == DiffusionReaction2d::ProblemA)
  {
    return RetType(impldiffusionreaction2d::TagProblemA{},
		   meshObj,
		   viscRecEn,
		   impldiffusionreaction2d::DefaultSourceProblemA<scalar_t>(),
		   0.01, 0.01);
  }

  else if (problemEnum == DiffusionReaction2d::GrayScott)
  {
    return RetType(impldiffusionreaction2d::TagProblemGrayScott{},
		   meshObj,
		   viscRecEn,
		   0.0002, 0.00005, 0.042, 0.062);
  }

  else{
    throw std::runtime_error("2D diffusion-reaction: invalid problem enum");
  }
}

template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<impldiffusionreaction2d::EigenApp<mesh_t>>
  >
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
RetType create_diffreac2d_problem_ov2(const mesh_t & meshObj,
				      DiffusionReaction2d problemEnum)
{
  return create_diffreac2d_problem_ov1(meshObj, problemEnum,
				       ViscousFluxReconstruction::FirstOrder);
}
#else
RetType create_problem_eigen(const mesh_t & meshObj,
			     DiffusionReaction2d problemEnum)
{
  return create_problem_eigen(meshObj, problemEnum, ViscousFluxReconstruction::FirstOrder);
}
#endif


// ----------------------------------------------------------
// create problem A with custom diffusion/reaction coeffs
// ----------------------------------------------------------
template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<impldiffusionreaction2d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
create_diffreac2d_problem_A_eigen_ov1
#else
create_diffusion_reaction_2d_problem_A_eigen
#endif
(const mesh_t & meshObj,
 ViscousFluxReconstruction viscFluxRecEnum,
 typename mesh_t::scalar_t diffusion,
 typename mesh_t::scalar_t reaction)
{
  using scalar_t = typename mesh_t::scalar_t;
  return RetType(impldiffusionreaction2d::TagProblemA{},
		 meshObj,
		 viscFluxRecEnum,
		 impldiffusionreaction2d::DefaultSourceProblemA<scalar_t>(),
		 diffusion, reaction);
}

// ----------------------------------------------------------
// create gray-scott with custom diffusion/reaction coeffs
// ----------------------------------------------------------
template<
  class mesh_t,
  class RetType = PublicProblemMixinCpp<impldiffusionreaction2d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
create_diffreac2d_problem_grayscott_eigen_ov1
#else
create_gray_scott_2d_problem_eigen
#endif
(const mesh_t & meshObj,
 ViscousFluxReconstruction viscFluxRecEnum,
 typename mesh_t::scalar_t diffusion_u,
 typename mesh_t::scalar_t diffusion_v,
 typename mesh_t::scalar_t feedRate,
 typename mesh_t::scalar_t killRate)
{
  using scalar_t = typename mesh_t::scalar_t;
  return RetType(impldiffusionreaction2d::TagProblemGrayScott{},
		 meshObj,
		 viscFluxRecEnum,
		 impldiffusionreaction2d::DefaultSourceProblemA<scalar_t>(),
		 diffusion_u, diffusion_v, feedRate, killRate);
}

} //end namespace pressiodemoapps
#endif


// template<class mesh_t>
// void check_stencil_admissibility(const mesh_t & meshObj,
// 			       ::pressiodemoapps::ViscousFluxReconstruction recEnum)
// {
//   const auto stencilSize = meshObj.stencilSize();
//   const auto check1 = stencilSizeCompatibleWithViscousFluxReconstruction(recEnum, stencilSize);
//   if (!check1){
//     throw std::runtime_error
//       ("Stencil size in the mesh object not compatible with desired viscous flux reconstruction.");
//   }
// }

// {
//   impldiffreac::check_stencil_admissibility(meshObj, recEnum);
//   return T(meshObj, probEnum, recEnum, diffusionCoeff, reactionCoeff);
// }

// template<class mesh_t, class T>
// T create_problem_for_pyC1d(const mesh_t & meshObj,
// 			   ::pressiodemoapps::DiffusionReaction1d probEnum,
// 			   ::pressiodemoapps::ViscousFluxReconstruction recEnum,
// 			   // source term is passed as a Python functor
// 			   pybind11::object pyFunctor,
// 			   typename mesh_t::scalar_t diffusionCoeff,
// 			   typename mesh_t::scalar_t reactionCoeff)
// {
//   impldiffreac::check_stencil_admissibility(meshObj, recEnum);
//   using scalar_t = typename mesh_t::scalar_t;
//   auto sourceWrapper = [=](const scalar_t & x, const scalar_t & evaltime, scalar_t & value){
//     value = pyFunctor.attr("__call__")(x, evaltime).template cast<scalar_t>();
//   };

//   return T(meshObj, probEnum, recEnum, sourceWrapper, diffusionCoeff, reactionCoeff);
// }

// template<class mesh_t, class T>
// T create_problem_for_pyC2d(const mesh_t & meshObj,
// 			   ::pressiodemoapps::DiffusionReaction2d probEnum,
// 			   ::pressiodemoapps::ViscousFluxReconstruction recEnum,
// 			   // source term is passed as a Python functor
// 			   pybind11::object pyFunctor,
// 			   typename mesh_t::scalar_t diffusionCoeff,
// 			   typename mesh_t::scalar_t reactionCoeff)
// {
//   impldiffreac::check_stencil_admissibility(meshObj, recEnum);
//   using scalar_t = typename mesh_t::scalar_t;
//   auto sourceWrapper = [=](const scalar_t & x,
// 			   const scalar_t & y,
// 			   const scalar_t & evaltime,
// 			   scalar_t & value){
//     value = pyFunctor.attr("__call__")(x, y, evaltime).template cast<scalar_t>();
//   };

//   return T(meshObj, probEnum, recEnum, sourceWrapper, diffusionCoeff, reactionCoeff);
// }

// template<class mesh_t, class T, class ProbType>
// T create_problem_for_pyD(const mesh_t & meshObj,
// 			 ProbType probEnum,
// 			 ::pressiodemoapps::ViscousFluxReconstruction recEnum,
// 			 typename mesh_t::scalar_t Da,
// 			 typename mesh_t::scalar_t Db,
// 			 typename mesh_t::scalar_t feedRate,
// 			 typename mesh_t::scalar_t killRate)
// {
//   impldiffreac::check_stencil_admissibility(meshObj, recEnum);
//   return T(meshObj, probEnum, recEnum, Da, Db, feedRate, killRate);
// }
//#endif
//} //end namespace impldiffreac
