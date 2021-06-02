
#ifndef PRESSIODEMOAPPS_ADVECTION_HPP_
#define PRESSIODEMOAPPS_ADVECTION_HPP_

#include <cmath>
#include <limits>

#include "mesh.hpp"
#include "weno.hpp"
#include "./impl/advection1d/fluxes.hpp"
#include "./impl/advection1d/linear_adv_impl.hpp"

namespace pressiodemoapps{

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class scalar_t, class mesh_t>
using PeriodicLinearAdvection1dEigen =
  pressiodemoapps::ad::impl::LinearAdvT<
  scalar_t, mesh_t,
  Eigen::Matrix<scalar_t,-1,1>,
  Eigen::Matrix<scalar_t,-1,1>,
  Eigen::Matrix<scalar_t,-1,1>
  >;
#endif

}
#endif
