
#ifndef PRESSIODEMOAPPS_MESH_HPP_
#define PRESSIODEMOAPPS_MESH_HPP_

#include <cmath>
#include <limits>
#include <fstream>
#include <array>
#include <iostream>
#include <sstream>

#include "./predicates/all.hpp"
#include "./container_fncs/all.hpp"
#include "./impl/mesh_ccu.hpp"

namespace pressiodemoapps{

template<class scalar_type = double>
auto load_cellcentered_uniform_mesh_eigen(const std::string & meshFilesPath)
{
  using return_type =  impl::CellCenteredUniformMesh<
    scalar_type,
    int32_t, // ordinal type
    Eigen::Matrix<scalar_type, -1, 1>, // coordinate storage type
    Eigen::Matrix<int32_t, -1, -1, Eigen::RowMajor> // graph type
    >;

  return return_type(meshFilesPath);
}

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class T>
T loadCellCenterUniformMesh(const std::string & meshFilesPath)
{
  return T(meshFilesPath);
}
#endif

}//end namespace pressiodemoapps
#endif
