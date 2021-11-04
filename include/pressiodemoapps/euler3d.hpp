
#ifndef PRESSIODEMOAPPS_EULER3D_INC_HPP_
#define PRESSIODEMOAPPS_EULER3D_INC_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./euler_compute_energy.hpp"

namespace pressiodemoapps{
enum class Euler3d{
  PeriodicSmooth,
  SedovSymmetry
};
}//end namespace pressiodemoapps

#include "./impl/euler_3d_prob_class.hpp"

namespace pressiodemoapps{
namespace impl{
template<class mesh_t, class T>
T createEe3dImpl(const mesh_t & meshObj,
		 ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		 ::pressiodemoapps::Euler3d probEnum,
		 ::pressiodemoapps::InviscidFluxScheme fluxEnum,
		 const int icId)
{
  const auto stencilSize = meshObj.stencilSize();
  const auto check1 = stencilSizeCompatibleWithInviscidFluxReconstruction(recEnum, stencilSize);
  if (!check1){
    throw std::runtime_error
      ("Stencil size in the mesh object not compatible with desired inviscid flux reconstruction.");
  }

  if (probEnum == ::pressiodemoapps::Euler3d::PeriodicSmooth)
  {
    if (!meshObj.isPeriodic()){
      throw std::runtime_error
      ("For euler3d::PeriodicSmooth, mesh must be periodic.");
    }
  }

  return T(meshObj, probEnum, recEnum, fluxEnum, icId);
}
} //end pressiodemoapps::impl

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler3d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum,
			::pressiodemoapps::InviscidFluxScheme fluxEnum,
			const int icId)
{
  using scalar_t = typename mesh_t::scalar_t;
  using p_t = ::pressiodemoapps::ee::impl::Euler3dAppT<
    scalar_t, mesh_t,
    Eigen::Matrix<scalar_t,-1,1>,
    Eigen::Matrix<scalar_t,-1,1>,
    Eigen::Matrix<scalar_t,-1,-1, Eigen::RowMajor>
    >;

  return impl::createEe3dImpl<mesh_t, p_t>(meshObj, recEnum, probEnum,
					   fluxEnum, icId);
}

template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Euler3d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return createProblemEigen(meshObj, probEnum, recEnum, InviscidFluxScheme::Rusanov, 1);
}
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createEuler3dForPyC(const mesh_t & meshObj,
		     ::pressiodemoapps::Euler3d probEnum,
		     ::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return impl::createEe3dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					 InviscidFluxScheme::Rusanov, 1);
}
#endif

}
#endif
