
#ifndef PRESSIODEMOAPPS_SWE2D_INC_HPP_
#define PRESSIODEMOAPPS_SWE2D_INC_HPP_

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
#include "Eigen/Core"
#endif

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./reconstruction.hpp"
#include "flux_enum.hpp"

namespace pressiodemoapps{
enum class Swe2d{
  GaussianPulse,
};
}//end namespace pressiodemoapps

#include "./impl_apps/swe2d/swe2d_app_impl.hpp"

namespace pressiodemoapps{
namespace impl{
template<class mesh_t, class T>
T createSwe2dImpl(const mesh_t & meshObj,
		 ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		 ::pressiodemoapps::Swe2d probEnum,
		 ::pressiodemoapps::InviscidFluxScheme fluxEnum)
{

  const auto stencilSize = meshObj.stencilSize();
  const auto check1 = stencilSizeCompatibleWithInviscidFluxReconstruction(recEnum, stencilSize);
  if (!check1){
    throw std::runtime_error
      ("Stencil size in the mesh object not compatible with desired inviscid flux reconstruction.");
  }

  return T(meshObj, probEnum, recEnum, fluxEnum);
}
} //end pressiodemoapps::impl

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Swe2d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum,
			::pressiodemoapps::InviscidFluxScheme fluxEnum)
{

  using scalar_t = typename mesh_t::scalar_t;
  using p_t = ::pressiodemoapps::swe::impl::Swe2dAppT<
    scalar_t, mesh_t,
    Eigen::Matrix<scalar_t,-1,1>,
    Eigen::Matrix<scalar_t,-1,1>,
    Eigen::Matrix<scalar_t,-1,-1, Eigen::RowMajor>
    >;

  return impl::createSwe2dImpl<mesh_t, p_t>(meshObj, recEnum, probEnum,
					   fluxEnum);
}

template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Swe2d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return createProblemEigen(meshObj, probEnum, recEnum,
			    InviscidFluxScheme::Rusanov);
}

#endif


#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createSwe2dForPyA(const mesh_t & meshObj,
		     ::pressiodemoapps::Swe2d probEnum,
		     ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		     ::pressiodemoapps::InviscidFluxScheme fluxEnum)
{
  return impl::createSwe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					 fluxEnum);
}


template<class mesh_t, class T>
T createSwe2dForPyB(const mesh_t & meshObj,
		     ::pressiodemoapps::Swe2d probEnum,
		     ::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return impl::createSwe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					 InviscidFluxScheme::Rusanov);
}


template<class mesh_t, class T>
T createSwe2dForPyC(const mesh_t & meshObj,
		     ::pressiodemoapps::Swe2d probEnum,
		     ::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return impl::createSwe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					 InviscidFluxScheme::Rusanov);
}

#endif


}
#endif
