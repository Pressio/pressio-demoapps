
#ifndef PRESSIODEMOAPPS_EULER1D_INC_HPP_
#define PRESSIODEMOAPPS_EULER1D_INC_HPP_

#include <cmath>
#include <limits>

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
#include "Eigen/Core"
#endif

#include "./predicates/all.hpp"
#include "./resize.hpp"
#include "./extent.hpp"
#include "./weno.hpp"
#include "./reconstruction_enums.hpp"
#include "./mesh.hpp"
#include "./impl/eulerCommon/energy.hpp"
#include "./impl/eulerCommon/fluxes.hpp"
#include "./impl/stencil_filler.hpp"
#include "./impl/reconstructor_from_stencil.hpp"
#include "./impl/reconstructor_from_state.hpp"

namespace pressiodemoapps{
enum class euler1dproblemsEnum{
  // for periodic, the user must provide a periodic domain and IC
  periodic,

  // sod taken from: https://www.mdpi.com/2227-7390/6/10/211/pdf
  // x \in [-0.5, 0.5]
  // initial condition is provided by the app class
  sod
};
}

#include "./impl/euler1d/ghost_filler.hpp"
#include "./impl/euler1d/initial_condition.hpp"
#include "./impl/euler1d/euler1d_app_impl.hpp"

namespace pressiodemoapps{

namespace impl{
template<class mesh_t, class T>
T createEe1dImpl(const mesh_t & meshObj,
		    pressiodemoapps::reconstructionEnum enIn,
		    pressiodemoapps::euler1dproblemsEnum probEn)
{
  meshObj.checkStencilSupportsOrder(enIn);

  if (probEn == pressiodemoapps::euler1dproblemsEnum::periodic){
    if (!meshObj.isPeriodic()){
      throw std::runtime_error
      ("For periodic euler1d, mesh must be periodic.");
    }
  }

  return T(meshObj, enIn, probEn);
}
} //end impl

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t, class scalar_t = double>
using Euler1dEigen =
  pressiodemoapps::ee::impl::Euler1dAppT<
  scalar_t, mesh_t,
  Eigen::Matrix<scalar_t,-1,1>,
  Eigen::Matrix<scalar_t,-1,1>,
  Eigen::Matrix<scalar_t,-1,1>
  >;

template<class mesh_t>
auto createEuler1dEigen(const mesh_t & meshObj,
			pressiodemoapps::reconstructionEnum enIn,
			pressiodemoapps::euler1dproblemsEnum probEn)
{
  using p_t = Euler1dEigen<mesh_t>;
  return impl::createEe1dImpl<mesh_t, p_t>(meshObj, enIn, probEn);
}
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t, class T>
T createEuler1dForPy(const mesh_t & meshObj,
			pressiodemoapps::reconstructionEnum enIn,
			pressiodemoapps::euler1dproblemsEnum probEn)
{
  return impl::createEe1dImpl<mesh_t, T>(meshObj, enIn, probEn);
}
#endif

}//end namespace pressiodemoapps
#endif
