
#ifndef PRESSIODEMOAPPS_ADVECTION_HPP_
#define PRESSIODEMOAPPS_ADVECTION_HPP_

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./mesh.hpp"
#include "./schemes_info.hpp"

namespace pressiodemoapps{
enum class Advection1d{
  PeriodicLinear
};
}//end namespace pressiodemoapps

#include "./impl/advection_1d_prob_class.hpp"

namespace pressiodemoapps{
namespace impl{
template<class mesh_t, class T>
T createLinAdv1dImpl(const mesh_t & meshObj,
		     ::pressiodemoapps::Advection1d probEnum,
		     InviscidFluxReconstruction recEnum)
{
  const auto stencilSize = meshObj.stencilSize();
  const auto check1 = stencilSizeCompatibleWithInviscidFluxReconstruction(recEnum, stencilSize);
  if (!check1){
    throw std::runtime_error
      ("advection: stencil size in the mesh object not compatible with desired inviscid flux reconstruction.");
  }

  // the mesh should be periodic
  if (!meshObj.isPeriodic()){
    throw std::runtime_error("For periodic linear adv1d, mesh must be periodic.");
  }

  return T(meshObj, probEnum, recEnum);
}
} //end impl

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t>
auto createProblemEigen(const mesh_t & meshObj,
			::pressiodemoapps::Advection1d probEnum,
			::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  using scalar_t = typename mesh_t::scalar_t;
  using p_type = ::pressiodemoapps::impladv::EigenAdvection1dAppWithJacobian<scalar_t, mesh_t>;
  return impl::createLinAdv1dImpl<mesh_t, p_type>(meshObj, probEnum, recEnum);
}
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createAdv1dForPy(const mesh_t & meshObj,
		   ::pressiodemoapps::Advection1d probEnum,
		   ::pressiodemoapps::InviscidFluxReconstruction recEnum)
{
  return impl::createLinAdv1dImpl<mesh_t, T>(meshObj, probEnum, recEnum);
}
#endif

}//end namespace pressiodemoapps
#endif
