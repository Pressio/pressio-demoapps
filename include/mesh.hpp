
#ifndef PRESSIODEMOAPPS_MESH_HPP_
#define PRESSIODEMOAPPS_MESH_HPP_

#include <cmath>
#include <limits>

#include "./impl/mesh/help_fncs.hpp"
#include "./impl/mesh/ccu_mesh.hpp"

namespace pressiodemoapps{

template<class ... Args>
using CellCenteredUniformMeshEigen =
  pressiodemoapps::impl::CellCenteredUniformMeshEigen<Args...>;

}

#endif
