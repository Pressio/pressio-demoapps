
#ifndef PRESSIODEMOAPPS_MESH_HPP_
#define PRESSIODEMOAPPS_MESH_HPP_

#include <cmath>
#include <limits>
#include <fstream>
#include <array>
#include <iostream>
#include <sstream>

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
#include "./meta/is_native_pybind_array.hpp"
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
#include "Eigen/Core"
#endif

#include "./impl/mesh/read_info.hpp"
#include "./impl/mesh/read_coords.hpp"
#include "./impl/mesh/read_connectivity.hpp"
#include "./impl/mesh/help_fncs.hpp"
#include "./impl/mesh/ccu_mesh.hpp"

namespace pressiodemoapps{

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class scalar_type>
using CellCenteredUniformMeshEigen =
  pressiodemoapps::impl::CellCenteredUniformMesh<
  scalar_type,
  int32_t,
  Eigen::Matrix<scalar_type, -1, 1>,
  Eigen::Matrix<int32_t, -1, -1, Eigen::RowMajor>,
  false
  >;
#endif

}
#endif
