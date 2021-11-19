
#ifndef PRESSIODEMOAPPS_SWE2D_INC_HPP_
#define PRESSIODEMOAPPS_SWE2D_INC_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"
#include "./adapter_mixins.hpp"

namespace pressiodemoapps{
enum class Swe2d{
  SlipWall,
};
}//end namespace pressiodemoapps

#include "./impl/swe_2d_prob_class.hpp"

namespace pressiodemoapps{
namespace implswe2d{

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

#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createProblemForPyA(const mesh_t & meshObj,
		      ::pressiodemoapps::Swe2d probEnum,
		      ::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return createSwe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
					  InviscidFluxScheme::Rusanov, 1);
}

template<class mesh_t, class T>
T createProblemForPyB(const mesh_t & meshObj,
		      ::pressiodemoapps::Swe2d probEnum,
		      ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		      const int ic)
{
  return createSwe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
				   InviscidFluxScheme::Rusanov, ic);
}
#endif

} //end pressiodemoapps::impl

template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Swe2d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum,
			::pressiodemoapps::InviscidFluxScheme fluxEnum,
			int icId)
{
  using p_t = ::pressiodemoapps::implswe::EigenSwe2dApp<mesh_t>;
  using RetType = PublicProblemMixinCpp<p_t>;
  return implswe2d::createSwe2dImpl<mesh_t, RetType>(meshObj, recEnum,
						     probEnum, fluxEnum, icId);
}

template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Swe2d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return createProblemEigen(meshObj, probEnum, recEnum,
				    InviscidFluxScheme::Rusanov, 1);
}

}//end namespace pressiodemoapps
#endif
