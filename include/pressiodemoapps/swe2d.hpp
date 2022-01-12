
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

template<class mesh_t, class T, class ...Args>
T createSwe2dImpl(const mesh_t & meshObj,
		  ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		  ::pressiodemoapps::Swe2d probEnum,
		  ::pressiodemoapps::InviscidFluxScheme fluxEnum,
		  Args && ... args)
{

  const auto stencilSize = meshObj.stencilSize();
  const auto check1 = stencilSizeCompatibleWithInviscidFluxReconstruction(recEnum, stencilSize);
  if (!check1){
    throw std::runtime_error
      ("Stencil size in the mesh object not compatible with desired inviscid flux reconstruction.");
  }

  return T(meshObj, probEnum, recEnum, fluxEnum, std::forward<Args>(args)...);
}

#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T create_problem_for_pyA(const mesh_t & meshObj,
			 ::pressiodemoapps::Swe2d probEnum,
			 ::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return createSwe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
				    InviscidFluxScheme::Rusanov, 1);
}

template<class mesh_t, class T>
T create_problem_for_pyB(const mesh_t & meshObj,
			 ::pressiodemoapps::Swe2d probEnum,
			 ::pressiodemoapps::InviscidFluxReconstruction recEnum,
			 typename mesh_t::scalar_t gravity,
			 typename mesh_t::scalar_t coriolis)
{
  return createSwe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
				    InviscidFluxScheme::Rusanov,
				    gravity, coriolis);
}

template<class mesh_t, class T>
T create_problem_for_pyC(const mesh_t & meshObj,
			 ::pressiodemoapps::Swe2d probEnum,
			 ::pressiodemoapps::InviscidFluxReconstruction recEnum,
			 typename mesh_t::scalar_t gravity,
			 typename mesh_t::scalar_t coriolis,
			 typename mesh_t::scalar_t pulseMagnitude)
{
  return createSwe2dImpl<mesh_t, T>(meshObj, recEnum, probEnum,
				    InviscidFluxScheme::Rusanov,
				    gravity, coriolis, pulseMagnitude);
}
#endif

} //end pressiodemoapps::impl

template<class mesh_t, class ...Args>
auto create_problem_eigen(const mesh_t & meshObj,
			  ::pressiodemoapps::Swe2d probEnum,
			  ::pressiodemoapps::InviscidFluxReconstruction recEnum,
			  Args && ... args)
{

  using p_t = ::pressiodemoapps::implswe::EigenSwe2dApp<mesh_t>;
  using RetType = PublicProblemMixinCpp<p_t>;
  return implswe2d::createSwe2dImpl<mesh_t, RetType>(meshObj, recEnum, probEnum,
						     InviscidFluxScheme::Rusanov,
						     std::forward<Args>(args)...);
}

}//end namespace pressiodemoapps
#endif
