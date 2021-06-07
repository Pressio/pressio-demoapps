
#ifndef PRESSIODEMOAPPS_MESH_HPP_
#define PRESSIODEMOAPPS_MESH_HPP_

#include <cmath>
#include <limits>
#include <fstream>
#include <array>
#include <iostream>
#include <sstream>

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
#include "Eigen/Core"
#endif

#include "./predicates/all.hpp"
#include "./resize.hpp"
#include "./extent.hpp"
#include "./reconstruction_enums.hpp"
#include "./impl/mesh/ccu_mesh.hpp"

namespace pressiodemoapps{

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class scalar_type = double>
using CellCenteredUniformMeshEigen =
  pressiodemoapps::impl::CellCenteredUniformMesh<
  scalar_type,
  int32_t,
  Eigen::Matrix<scalar_type, -1, 1>,
  Eigen::Matrix<int32_t, -1, -1, Eigen::RowMajor>,
  false
  >;

auto loadCellCenterUniformMeshEigen(const std::string & meshFilesPath)
{
  return CellCenteredUniformMeshEigen<double>(meshFilesPath);
}
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class T>
T loadCellCenterUniformMesh(const std::string & meshFilesPath)
{
  return T(meshFilesPath);
}
#endif

}
#endif
