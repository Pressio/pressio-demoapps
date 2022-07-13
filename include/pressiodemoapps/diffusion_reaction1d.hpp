
#ifndef PRESSIODEMOAPPS_DIFFUSION_REACTION_1D_HPP_
#define PRESSIODEMOAPPS_DIFFUSION_REACTION_1D_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./adapter_cpp.hpp"
#include "./adapter_py.hpp"
#include "mypi.hpp"

namespace pressiodemoapps{

// ----------------------------------------------------------
// enums identifying the problems
// ----------------------------------------------------------
enum class DiffusionReaction1d
{
  ProblemA
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
};

// ----------------------------------------------------------
// default source terms
// ----------------------------------------------------------
namespace impldiffusionreaction1d{
template<class scalar_t>
struct DefaultSourceProblemA{
  void operator()(const scalar_t & x,
		  const scalar_t & evaltime,
		  scalar_t & value)
  {
    (void) evaltime;
    value = std::sin(M_PI*x) * x*x * 4.*std::cos(4.*M_PI*x);
  }
};
}//end namespace impldiffusionreaction1d
}//end namespace pressiodemoapps


// this include is here because needs visiblity of the enums above
#include "./impl/diffusion_reaction_1d_prob_class.hpp"


namespace pressiodemoapps{

// ----------------------------------------------------------
// create default problem
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impldiffusionreaction1d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// bindings need unique naming or we get error associated with overloads
create_diffusion_reaction_1d_problem_default_for_py
#else
create_problem_eigen
#endif
(const mesh_t & meshObj,
 DiffusionReaction1d problemEnum)
{

  using scalar_t = typename mesh_t::scalar_t;
  if (problemEnum == DiffusionReaction1d::ProblemA)
  {
    return RetType(impldiffusionreaction1d::TagProblemA{},
		   meshObj,
		   ViscousFluxReconstruction::FirstOrder,
		   impldiffusionreaction1d::DefaultSourceProblemA<scalar_t>(),
		   0.01, 0.01);
  }
  else{
    throw std::runtime_error("1D diffusion-reaction: invalid problem enum");
  }
}

// ----------------------------------------------------------
// create problem A with custom diffusion/reaction coeffs
// ----------------------------------------------------------

template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impldiffusionreaction1d::EigenApp<mesh_t>>
  >
RetType
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
// bindings need unique naming or we get error associated with overloads
create_diffusion_reaction_1d_problem_A_ov1_for_py
#else
create_diffusion_reaction_1d_problem_A_eigen
#endif
(const mesh_t & meshObj,
 typename mesh_t::scalar_t diffusion,
 typename mesh_t::scalar_t reaction)
{

  using scalar_t = typename mesh_t::scalar_t;
  return RetType(impldiffusionreaction1d::TagProblemA{},
		 meshObj,
		 ViscousFluxReconstruction::FirstOrder,
		 impldiffusionreaction1d::DefaultSourceProblemA<scalar_t>(),
		 diffusion, reaction);
}

// ----------------------------------------------------------
// create problem A with custom diffusion/reaction and source
// ----------------------------------------------------------
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<
  class mesh_t,
  class RetType = PublicProblemEigenMixinCpp<impldiffusionreaction1d::EigenApp<mesh_t>>
  >
RetType create_diffusion_reaction_1d_problem_A_ov2_for_py(const mesh_t & meshObj,
							  // source term is passed as a Python functor
							  pybind11::object pyFunctor,
							  typename mesh_t::scalar_t diffusion,
							  typename mesh_t::scalar_t reaction)
{
  using scalar_t = typename mesh_t::scalar_t;
  auto sourceWrapper = [=](const scalar_t & x,
			   const scalar_t & evaltime,
			   scalar_t & value){
    value = pyFunctor.attr("__call__")(x, evaltime).template cast<scalar_t>();
  };

  return RetType(impldiffusionreaction1d::TagProblemA{},
		 meshObj, ViscousFluxReconstruction::FirstOrder,
		 sourceWrapper, diffusion, reaction);
}

#else

template<
  class mesh_t,
  class SourceType,
  class RetType = PublicProblemEigenMixinCpp<impldiffusionreaction1d::EigenApp<mesh_t>>
  >
RetType create_diffusion_reaction_1d_problem_A_eigen(const mesh_t & meshObj,
						     SourceType sourceF,
						     typename mesh_t::scalar_t diffusion,
						     typename mesh_t::scalar_t reaction)
{
  return RetType(impldiffusionreaction1d::TagProblemA{},
		 meshObj, ViscousFluxReconstruction::FirstOrder,
		 sourceF, diffusion, reaction);
}
#endif

} //end namespace pressiodemoapps
#endif
