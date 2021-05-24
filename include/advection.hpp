
#ifndef PRESSIODEMOAPPS_ADVECTION_HPP_
#define PRESSIODEMOAPPS_ADVECTION_HPP_

#include <cmath>
#include <limits>

#include "mesh.hpp"
#include "weno.hpp"
#include "./impl/advection1d/fluxes.hpp"

#include "Eigen/Core"
#include "Eigen/Sparse"
#include "./impl/advection1d/eigen_linear.hpp"

namespace pressiodemoapps{

template<class scalar_t, class mesh_t>
using PeriodicLinearAdvection1dEigen =
  pressiodemoapps::ad::impl::EigenLinearAdv<scalar_t, mesh_t>;

}

#endif
