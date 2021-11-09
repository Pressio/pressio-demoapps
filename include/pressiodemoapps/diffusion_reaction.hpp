
#ifndef PRESSIODEMOAPPS_DIFFUSION_REACTION_HPP_
#define PRESSIODEMOAPPS_DIFFUSION_REACTION_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"

namespace pressiodemoapps{
enum class DiffusionReaction1d{
  /*
    ds/dt = D d^2s/d^2x + k*s^2 + u(x, t)

    Default:
    D,k = 0.01, 0.01
    u(s) = sin(M_PI*x) * x*x * 4.*std::cos(4.*M_PI*x);

    D,k,u(s) can be provided to constructor
    u(s) must be a functor:
      void operator()(const scalar_t & x, const scalar_t & time, scalar_t & value);

    BC: ghost cells are set such that s is zero at boundary
   */
  ProblemA
};
}//end namespace pressiodemoapps

#include "./impl/diffusion_reaction_1d_prob_class.hpp"

namespace pressiodemoapps{
namespace impldiffreac{

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createDiffReac1dForPyA(const mesh_t & meshObj,
			 ::pressiodemoapps::DiffusionReaction1d probEnum,
			 ::pressiodemoapps::ViscousFluxReconstruction recEnum)
{
  return T(meshObj, probEnum, recEnum);
}

template<class mesh_t, class T>
T createDiffReac1dForPyB(const mesh_t & meshObj,
			 ::pressiodemoapps::DiffusionReaction1d probEnum,
			 ::pressiodemoapps::ViscousFluxReconstruction recEnum,
			 typename mesh_t::scalar_t diffusionCoeff,
			 typename mesh_t::scalar_t reactionCoeff)
{
  return T(meshObj, probEnum, recEnum, diffusionCoeff, reactionCoeff);
}

template<class mesh_t, class T>
T createDiffReac1dForPyC(const mesh_t & meshObj,
			 ::pressiodemoapps::DiffusionReaction1d probEnum,
			 ::pressiodemoapps::ViscousFluxReconstruction recEnum,
			 pybind11::object pyFunctor,
			 typename mesh_t::scalar_t diffusionCoeff,
			 typename mesh_t::scalar_t reactionCoeff)
{
  using scalar_t = typename mesh_t::scalar_t;
  auto sourceWrapper = [=](const scalar_t & x, const scalar_t & evaltime, scalar_t & value){
    value = pyFunctor.attr("__call__")(x, evaltime).template cast<scalar_t>();
  };

  return T(meshObj, probEnum, recEnum, sourceWrapper, diffusionCoeff, reactionCoeff);
}
#endif
}//end namespace impldiffreac


#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t, class ...Args>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::DiffusionReaction1d probEnum,
			::pressiodemoapps::ViscousFluxReconstruction recEnum,
			Args && ... args)
{
  const auto stencilSize = meshObj.stencilSize();
  const auto check1 = stencilSizeCompatibleWithViscousFluxReconstruction(recEnum, stencilSize);
  if (!check1){
    throw std::runtime_error
      ("Stencil size in the mesh object not compatible with desired viscous flux reconstruction.");
  }

  using scalar_t = typename mesh_t::scalar_t;
  using return_type = ::pressiodemoapps::impldiffreac::DiffReac1dApp<
    scalar_t, mesh_t,
    Eigen::Matrix<scalar_t,-1,1>, // state
    Eigen::Matrix<scalar_t,-1,1>, // velocity
    Eigen::Matrix<scalar_t,-1,1>,  // ghost (in 1d, these are just vectors)
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor, int>
    >;
  return return_type(meshObj, probEnum, recEnum, std::forward<Args>(args)...);
}
#endif


}//end namespace pressiodemoapps
#endif
