/*
//@HEADER
// ************************************************************************
//
// mesh_ccu.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef PRESSIODEMOAPPS_CCU_MESH_HPP_
#define PRESSIODEMOAPPS_CCU_MESH_HPP_

#include "mesh_read_info.hpp"
#include "mesh_read_coords.hpp"
#include "mesh_read_connectivity.hpp"
#include <tuple>
#include <vector>

namespace pressiodemoapps{ namespace impl{

template<
  class ScalarType,
  class IndexType,
  class xyz_type,
  class mesh_g_type,
  bool is_binding = false
  >
class CellCenteredUniformMesh
{

public:
  // these are here for backward comp but need to remove
  using scalar_t  = ScalarType;
  using index_t	  = IndexType;
  using graph_t   = mesh_g_type;

  using scalar_type  = ScalarType;
  using index_type  = IndexType;
  using x_coords_type = xyz_type;
  using y_coords_type = xyz_type;
  using z_coords_type = xyz_type;

private:
  using indices_v_t = std::vector<index_type>;

public:
  CellCenteredUniformMesh() = default;

#if not defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  template<
    bool _is_binding = is_binding,
    typename std::enable_if<!_is_binding>::type * = nullptr
    >
  explicit CellCenteredUniformMesh(const std::string & meshDir)
  {
    allocateAndSetup(meshDir);
    m_meshIsFullyPeriodic = checkIfFullyPeriodic();
  }

#else
  // note that when doing bindings, I need to first construct
  // with {1,1} just so that the numpy array picks up these
  // are 2dim array otherwise it thinks they are 1d array.
  template<
    bool _is_binding = is_binding,
    typename std::enable_if<_is_binding>::type * = nullptr
    >
  explicit CellCenteredUniformMesh(const std::string & meshDir)
    : m_x(1), m_y(1), m_z(1), m_graph({1,1})
  {
    allocateAndSetup(meshDir);
    m_meshIsFullyPeriodic = checkIfFullyPeriodic();
  }
#endif

  int dimensionality() const{ return m_dim; }
  index_type stencilMeshSize() const{ return m_stencilMeshSize; }
  index_type sampleMeshSize() const{ return m_sampleMeshSize; }
  int stencilSize() const { return m_stencilSize; }

  const graph_t & graph() const{ return m_graph; }

  ScalarType dx() const{ return m_cellDeltas[0]; }
  ScalarType dy() const{ return m_cellDeltas[1]; }
  ScalarType dz() const{ return m_cellDeltas[2]; }
  ScalarType dxInv() const{ return m_cellDeltasInv[0]; }
  ScalarType dyInv() const{ return m_cellDeltasInv[1]; }
  ScalarType dzInv() const{ return m_cellDeltasInv[2]; }

  const x_coords_type & viewX() const{ return m_x; }
  const y_coords_type & viewY() const{ return m_y; }
  const z_coords_type & viewZ() const{ return m_z; }

  // number of cells that do NOT get close to any boundary
  // i.e. those for which the farthest stencil cell is still within the BD
  auto numCellsInner() const{ return m_rowsForCellsInner.size(); }

  // the rows of the graph that pertain "inner" cells
  const indices_v_t & graphRowsOfCellsAwayFromBd() const{
    return m_rowsForCellsInner;
  }

  // all the other cells, i.e. those that are close enough to a BD
  auto numCellsBd() const   { return m_rowsForCellsBd.size(); }

  // the rows of the graph that pertain "non-inner" cells
  const indices_v_t & graphRowsOfCellsNearBd() const{
    return m_rowsForCellsBd;
  }

  bool isFullyPeriodic() const{
    return m_meshIsFullyPeriodic;
  }

  // 1d
  bool hasBdLeft1d(const index_type rowInd) const{
    if (m_graph(rowInd, 1)==-1) return true;

    if (m_stencilSize>=5){
      if (m_graph(rowInd, 3)==-1) return true;
    }
    if (m_stencilSize==7){
      if (m_graph(rowInd, 5)==-1) return true;
    }
    return false;
  }

  bool hasBdRight1d(const index_type rowInd) const{
    if (m_graph(rowInd, 2)==-1) return true;

    if (m_stencilSize>=5){
      if (m_graph(rowInd, 4)==-1) return true;
    }
    if (m_stencilSize==7){
      if (m_graph(rowInd, 6)==-1) return true;
    }
    return false;
  }

  // 2d
  bool hasBdLeft2d(const index_type rowInd) const{
    if (m_graph(rowInd, 1)==-1) {
      return true;
    }

    if (m_stencilSize>=5){
      if (m_graph(rowInd, 5)==-1) return true;
    }
    if (m_stencilSize==7){
      if (m_graph(rowInd, 9)==-1) return true;
    }
    return false;
  }

  bool hasBdFront2d(const index_type rowInd) const{
    if (m_graph(rowInd, 2)==-1) {
      return true;
    }

    if (m_stencilSize>=5){
      if (m_graph(rowInd, 6)==-1) return true;
    }
    if (m_stencilSize==7){
      if (m_graph(rowInd, 10)==-1) return true;
    }
    return false;
  }

  bool hasBdRight2d(const index_type rowInd) const{
    if (m_graph(rowInd, 3)==-1) {
      return true;
    }

    if (m_stencilSize>=5){
      if (m_graph(rowInd, 7)==-1) return true;
    }
    if (m_stencilSize==7){
      if (m_graph(rowInd, 11)==-1) return true;
    }
    return false;
  }

  bool hasBdBack2d(const index_type rowInd) const{
    if (m_graph(rowInd, 4)==-1) {
      return true;
    }

    if (m_stencilSize>=5){
      if (m_graph(rowInd, 8)==-1) return true;
    }
    if (m_stencilSize==7){
      if (m_graph(rowInd, 12)==-1) return true;
    }
    return false;
  }

  // 3d
  bool hasBdLeft3d(const index_type rowInd) const{
    if (m_graph(rowInd, 1)==-1) return true;

    if (m_stencilSize==5){
      if (m_graph(rowInd, 7)==-1) return true;
    }
    return false;
  }

  bool hasBdFront3d(const index_type rowInd) const{
    if (m_graph(rowInd, 2)==-1) return true;

    if (m_stencilSize==5){
      if (m_graph(rowInd, 8)==-1) return true;
    }
    return false;
  }

  bool hasBdRight3d(const index_type rowInd) const{
    if (m_graph(rowInd, 3)==-1) return true;

    if (m_stencilSize==5){
      if (m_graph(rowInd, 9)==-1) return true;
    }
    return false;
  }

  bool hasBdBack3d(const index_type rowInd) const{
    if (m_graph(rowInd, 4)==-1) return true;

    if (m_stencilSize==5){
      if (m_graph(rowInd, 10)==-1) return true;
    }
    return false;
  }

  bool hasBdBottom3d(const index_type rowInd) const{
    if (m_graph(rowInd, 5)==-1) return true;

    if (m_stencilSize==5){
      if (m_graph(rowInd, 11)==-1) return true;
    }
    return false;
  }

  bool hasBdTop3d(const index_type rowInd) const{
    if (m_graph(rowInd, 6)==-1) return true;

    if (m_stencilSize==5){
      if (m_graph(rowInd, 12)==-1) return true;
    }
    return false;
  }


private:
  // mesh connectivity is periodic if for each mesh cell,
  // all neighbors global GIDs are non-negative.
  bool checkIfFullyPeriodic(){
    // Note that we need to check along ALL directions.
    const auto n = (m_stencilSize-1)*m_dim;
    for (index_type i=0; i<m_sampleMeshSize; ++i){
      for (index_type j=0; j<n+1; ++j){
	if (m_graph(i, j) < 0){
	  return false;
	}
      }
    }
    return true;
  }

  auto boundsImpl(int i) const
  {
    using numlimits = std::numeric_limits<ScalarType>;
    auto res = std::make_tuple(numlimits::max(), numlimits::min());

    const auto & a = (i==0) ? m_x : (i==1) ? m_y : m_z;
    for (int i=0; i<m_stencilMeshSize; ++i){
      auto & v1 = res.template get<0>(res);
      auto & v2 = res.template get<1>(res);
      v1 = std::min(v1, a(i));
      v2 = std::max(v2, a(i));
    }
    return res;
  }

  void allocateAndSetup(const std::string & meshDir)
  {
    pressiodemoapps::impl::read_mesh_info(meshDir, m_dim,
					  m_cellDeltas,
					  m_cellDeltasInv,
					  m_stencilSize,
					  m_stencilMeshSize,
					  m_sampleMeshSize);
    if(m_dim==3 and m_stencilSize==7){
      throw std::runtime_error("3D mesh with 7pt stencil not yet supported.");
    }

    pressiodemoapps::resize(m_x, m_stencilMeshSize);
    pressiodemoapps::resize(m_y, m_stencilMeshSize);
    pressiodemoapps::resize(m_z, m_stencilMeshSize);
    pressiodemoapps::impl::read_mesh_coordinates(meshDir, m_dim, m_x, m_y, m_z);

    // compute the number of of neighbors
    const auto numNeighbors = (m_stencilSize-1)*m_dim;
    // the graph # of cols is num of neighbors + 1 becuase
    // the graph col(0) contains the global ID of the cell itself
    const auto graphNumCols = numNeighbors + 1;
    pressiodemoapps::resize(m_graph, m_sampleMeshSize, graphNumCols);
    pressiodemoapps::impl::read_mesh_connectivity(meshDir, m_graph, graphNumCols);

    // figure out how many cells are near the boundaries
    for (index_type it=0; it<m_sampleMeshSize; ++it)
    {
      if (m_dim==1)
      {
	const auto b1 = hasBdLeft1d(it);
	const auto b2 = hasBdRight1d(it);
	// if either is true, this cell is near the BD, so
	// add its graph row to corresponding list
	if (b1 or b2){
	  m_rowsForCellsBd.push_back(it);
	}
	else{
	  m_rowsForCellsInner.push_back(it);
	}
      }

      else if (m_dim==2)
      {
	const auto b1 = hasBdLeft2d(it);
	const auto b2 = hasBdFront2d(it);
	const auto b3 = hasBdRight2d(it);
	const auto b4 = hasBdBack2d(it);

	// if either is true, this cell is near the BD, so
	// add its graph row to corresponding list
	if (b1 or b2 or b3 or b4){
	  m_rowsForCellsBd.push_back(it);
	}
	else{
	  m_rowsForCellsInner.push_back(it);
	}
      }

      else if (m_dim==3)
      {
	const auto b1 = hasBdLeft3d(it);
	const auto b2 = hasBdFront3d(it);
	const auto b3 = hasBdRight3d(it);
	const auto b4 = hasBdBack3d(it);
	const auto b5 = hasBdBottom3d(it);
	const auto b6 = hasBdTop3d(it);

	// if either is true, this cell is near the BD, so
	// add its graph row to corresponding list
	if (b1 or b2 or b3 or b4 or b5 or b6){
	  m_rowsForCellsBd.push_back(it);
	}
	else{
	  m_rowsForCellsInner.push_back(it);
	}
      }
      else{
	throw std::runtime_error("Invalid dimension");
      }
    }
  }

private:
  int m_dim = {};
  std::array<ScalarType,3> m_cellDeltas{};
  std::array<ScalarType,3> m_cellDeltasInv{};
  int m_stencilSize = {};
  index_type m_stencilMeshSize = {};
  index_type m_sampleMeshSize  = {};
  x_coords_type m_x = {};
  y_coords_type m_y = {};
  z_coords_type m_z = {};

  /*
    graph:

    first col : contains GIDs of cells where we want velocity
    col 1,... : contains GIDs of neighboring cells needed for stencil
  */
  graph_t m_graph = {};

  indices_v_t m_rowsForCellsInner;
  indices_v_t m_rowsForCellsBd;

  bool m_meshIsFullyPeriodic = false;
};

}}
#endif
