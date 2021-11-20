
#ifndef PRESSIODEMOAPPS_CCU_MESH_HPP_
#define PRESSIODEMOAPPS_CCU_MESH_HPP_

#include "mesh_read_info.hpp"
#include "mesh_read_coords.hpp"
#include "mesh_read_connectivity.hpp"
#include <tuple>
#include <vector>

namespace pressiodemoapps{ namespace impl{

template<
  class scalar_type,
  class index_type,
  class xyz_type,
  class mesh_g_type,
  bool is_binding = false
  >
class CellCenteredUniformMesh
{

public:
  using scalar_t  = scalar_type;
  using index_t	  = index_type;
  using x_t	  = xyz_type;
  using y_t	  = xyz_type;
  using z_t	  = xyz_type;
  using graph_t   = mesh_g_type;
  using indices_v_t = std::vector<index_t>;

  CellCenteredUniformMesh() = delete;

#if not defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  template<
    bool _is_binding = is_binding,
    typename std::enable_if<!_is_binding>::type * = nullptr
    >
  explicit CellCenteredUniformMesh(const std::string & meshDir)
  {
    allocateAndSetup(meshDir);
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
  }
#endif

  int dimensionality() const{ return m_dim; }
  index_t stencilMeshSize() const{ return m_stencilMeshSize; }
  index_t sampleMeshSize() const{ return m_sampleMeshSize; }
  const int stencilSize() const { return m_stencilSize; }

  const graph_t & graph() const{ return m_graph; }

  auto boundsX() const{ return boundsImpl(0); }
  auto boundsY() const{ return boundsImpl(1); }
  auto boundsZ() const{ return boundsImpl(2); }

  const scalar_type dx() const{ return m_cellDeltas[0]; }
  const scalar_type dy() const{ return m_cellDeltas[1]; }
  const scalar_type dz() const{ return m_cellDeltas[2]; }
  const scalar_type dxInv() const{ return m_cellDeltasInv[0]; }
  const scalar_type dyInv() const{ return m_cellDeltasInv[1]; }
  const scalar_type dzInv() const{ return m_cellDeltasInv[2]; }

  const x_t & viewX() const{ return m_x; }
  const y_t & viewY() const{ return m_y; }
  const z_t & viewZ() const{ return m_z; }

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

  // mesh connectivity is periodic if for each mesh cell,
  // all neighbors global GIDs are non-negative.
  bool isPeriodic() const{
    // Note that we need to check along ALL directions.
    const auto n = (m_stencilSize-1)*m_dim;

    for (index_t i=0; i<m_sampleMeshSize; ++i){
      for (index_t j=0; j<n+1; ++j){
	if (m_graph(i, j) < 0){
	  return false;
	}
      }
    }
    return true;
  }

  // 1d
  bool hasBdLeft1d(const index_t rowInd) const{
    if (m_graph(rowInd, 1)==-1) return true;

    if (m_stencilSize>=5){
      if (m_graph(rowInd, 3)==-1) return true;
    }
    if (m_stencilSize==7){
      if (m_graph(rowInd, 5)==-1) return true;
    }
    return false;
  }

  bool hasBdRight1d(const index_t rowInd) const{
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
  bool hasBdLeft2d(const index_t rowInd) const{
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

  bool hasBdFront2d(const index_t rowInd) const{
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

  bool hasBdRight2d(const index_t rowInd) const{
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

  bool hasBdBack2d(const index_t rowInd) const{
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
  bool hasBdLeft3d(const index_t rowInd) const{
    if (m_graph(rowInd, 1)==-1) return true;

    if (m_stencilSize==5){
      if (m_graph(rowInd, 7)==-1) return true;
    }
    return false;
  }

  bool hasBdFront3d(const index_t rowInd) const{
    if (m_graph(rowInd, 2)==-1) return true;

    if (m_stencilSize==5){
      if (m_graph(rowInd, 8)==-1) return true;
    }
    return false;
  }

  bool hasBdRight3d(const index_t rowInd) const{
    if (m_graph(rowInd, 3)==-1) return true;

    if (m_stencilSize==5){
      if (m_graph(rowInd, 9)==-1) return true;
    }
    return false;
  }

  bool hasBdBack3d(const index_t rowInd) const{
    if (m_graph(rowInd, 4)==-1) return true;

    if (m_stencilSize==5){
      if (m_graph(rowInd, 10)==-1) return true;
    }
    return false;
  }

  bool hasBdBottom3d(const index_t rowInd) const{
    if (m_graph(rowInd, 5)==-1) return true;

    if (m_stencilSize==5){
      if (m_graph(rowInd, 11)==-1) return true;
    }
    return false;
  }

  bool hasBdTop3d(const index_t rowInd) const{
    if (m_graph(rowInd, 6)==-1) return true;

    if (m_stencilSize==5){
      if (m_graph(rowInd, 12)==-1) return true;
    }
    return false;
  }


private:
  auto boundsImpl(int i) const
  {
    using numlimits = std::numeric_limits<scalar_type>;
    auto res = std::make_tuple(numlimits::max(), numlimits::min());

    const auto & a = (i==0) ? m_x : (i==1) ? m_y : m_z;
    for (int i=0; i<m_stencilMeshSize; ++i){
      auto & v1 = res.get<0>(res);
      auto & v2 = res.get<1>(res);
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
    for (index_t it=0; it<m_sampleMeshSize; ++it)
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
  std::array<scalar_type,3> m_cellDeltas{};
  std::array<scalar_type,3> m_cellDeltasInv{};
  int m_stencilSize = {};
  index_t m_stencilMeshSize = {};
  index_t m_sampleMeshSize  = {};
  x_t m_x = {};
  y_t m_y = {};
  z_t m_z = {};

  /*
    graph:

    first col : contains GIDs of cells where we want velocity
    col 1,... : contains GIDs of neighboring cells needed for stencil
  */
  graph_t m_graph = {};

  indices_v_t m_rowsForCellsInner;
  indices_v_t m_rowsForCellsBd;

};

}}
#endif
