
#ifndef PRESSIODEMOAPPS_SWE2D_INC_HPP_
#define PRESSIODEMOAPPS_SWE2D_INC_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"

namespace pressiodemoapps{
enum class Swe2d{
  SlipWall,
};
}//end namespace pressiodemoapps

#include "./impl/swe_2d_prob_class.hpp"

namespace pressiodemoapps{
namespace impl{
template<class mesh_t, class T>
T createSwe2dImpl(const mesh_t & meshObj,
		 ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		 ::pressiodemoapps::Swe2d probEnum,
		 ::pressiodemoapps::InviscidFluxScheme fluxEnum,
		  const int icId)
{

  const auto stencilSize = meshObj.stencilSize();
  const auto check1 = stencilSizeCompatibleWithInviscidFluxReconstruction(recEnum, stencilSize);
  if (!check1){
    throw std::runtime_error
      ("Stencil size in the mesh object not compatible with desired inviscid flux reconstruction.");
  }

  return T(meshObj, probEnum, recEnum, fluxEnum,icId);
}
} //end pressiodemoapps::impl

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Swe2d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum,
			::pressiodemoapps::InviscidFluxScheme fluxEnum,
      int initCondIdentifier)
{

  using scalar_t = typename mesh_t::scalar_t;
  using p_t = ::pressiodemoapps::implswe::Swe2dAppT<
    scalar_t, mesh_t,
    Eigen::Matrix<scalar_t,-1,1>,
    Eigen::Matrix<scalar_t,-1,1>,
    Eigen::Matrix<scalar_t,-1,-1, Eigen::RowMajor>,
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor, int>
    >;

  return impl::createSwe2dImpl<mesh_t, p_t>(meshObj, recEnum, probEnum,
					   fluxEnum,initCondIdentifier);
}

template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Swe2d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return createProblemEigen(meshObj, probEnum, recEnum,
			    InviscidFluxScheme::Rusanov,1);
}
#endif


#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createSwe2dForPyA(const mesh_t & meshObj,
		     ::pressiodemoapps::Swe2d probEnum,
		     ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		     ::pressiodemoapps::InviscidFluxScheme fluxEnum,
         const int initCondIdentifier = 1)
{
  return impl::createSwe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					 fluxEnum,initCondIdentifier);
}


template<class mesh_t, class T>
T createSwe2dForPyB(const mesh_t & meshObj,
		     ::pressiodemoapps::Swe2d probEnum,
		     ::pressiodemoapps::InviscidFluxReconstruction recEnum,
         const int icId = 1)
{
  return impl::createSwe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					 InviscidFluxScheme::Rusanov,icId);
}


template<class mesh_t, class T>
T createSwe2dForPyC(const mesh_t & meshObj,
		     ::pressiodemoapps::Swe2d probEnum,
		     ::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return impl::createSwe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					 InviscidFluxScheme::Rusanov,1);
}
#endif

}
#endif
