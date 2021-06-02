
#ifndef PRESSIODEMOAPPS_EULER_HPP_
#define PRESSIODEMOAPPS_EULER_HPP_

#include <cmath>
#include <limits>

#include "mesh.hpp"
#include "weno.hpp"
#include "./impl/eulerCommon/energy.hpp"
#include "./impl/eulerCommon/fluxes.hpp"

#include "./impl/euler2d/initial_condition.hpp"

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "./impl/euler1d/eigen_app.hpp"
#include "./impl/euler2d/eigen_app.hpp"
#endif

namespace pressiodemoapps{

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class scalar_t, class mesh_t>
using PeriodicEuler1dEigen =
  pressiodemoapps::ee::impl::EigenApp1d<scalar_t, mesh_t, 0>;

template<class scalar_t, class mesh_t>
using Sod1dEigen =
  pressiodemoapps::ee::impl::EigenApp1d<scalar_t, mesh_t, 1>;

template<class scalar_t, class mesh_t>
using PeriodicEuler2dEigen =
  pressiodemoapps::ee::impl::EigenApp2d<scalar_t, mesh_t, 0>;

template<class scalar_t, class mesh_t>
using Sedov2dEigen =
  pressiodemoapps::ee::impl::EigenApp2d<scalar_t, mesh_t, 1>;

template<class scalar_t, class mesh_t>
using Riemann2dEigen =
  pressiodemoapps::ee::impl::EigenApp2d<scalar_t, mesh_t, 2>;
#endif

}

#endif
