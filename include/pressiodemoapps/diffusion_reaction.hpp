
#ifndef PRESSIODEMOAPPS_DIFFUSION_REACTION_HPP_
#define PRESSIODEMOAPPS_DIFFUSION_REACTION_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./flux_enum.hpp"
#include "./reconstruction.hpp"

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
   */
  ProblemA
};
}//end namespace pressiodemoapps

#include "./impl_apps/diffusionReaction1d/diff_reac_1d_impl.hpp"

namespace pressiodemoapps{
// namespace impldiffreac{
// template<class mesh_t, class T>
// T createDiffReac1d(const mesh_t & meshObj,
// 		     ::pressiodemoapps::DiffusionReaction1d probEnum)
// {
//   return T(meshObj, probEnum);
// }
// } //end pressiodemoapps::impldiffreac

// #ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
// template<class mesh_t, class T>
// T createAdv1dForPy(const mesh_t & meshObj,
// 		   ::pressiodemoapps::DiffusionReaction1d probEnum)
// {
//   return T(meshObj, probEnum);
// }
// #endif

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
  using return_type = ::pressiodemoapps::impldiffreac::DiffReac1dAppT<
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
