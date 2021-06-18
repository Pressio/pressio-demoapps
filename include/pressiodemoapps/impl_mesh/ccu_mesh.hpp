
#ifndef PRESSIODEMOAPPS_CCU_MESH_HPP_
#define PRESSIODEMOAPPS_CCU_MESH_HPP_

#include "./read_info.hpp"
#include "./read_coords.hpp"
#include "./read_connectivity.hpp"
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
  // with {1,1} just so that the numpy array picks up they
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

  // void checkStencilSupportsOrder(pressiodemoapps::reconstructionType eIn) const
  // {
  //   if (m_stencilSize == 3)
  //   {
  //     if (eIn == reconstructionType::Weno5){
  // 	throw std::runtime_error
  // 	  ("Mesh stencil size not large enough for target reconstruction order.");
  //     }
  //   }
  // }

  bool isPeriodic() const{
    const auto gSize = (m_stencilSize-1)*m_dim;

    // mesh connectivity is periodic if all neighboring GIDs are positive
    for (index_t i=0; i<m_sampleMeshSize; ++i){
      for (index_t j=0; j<gSize+1; ++j){
	if (m_graph(i, j) < 0){
	  return false;
	}
      }
    }
    return true;
  }

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

  auto numCellsInner() const{ return m_rowsForCellsInner.size(); }
  auto numCellsBd() const   { return m_rowsForCellsBd.size(); }

  const indices_v_t & graphRowsOfCellsAwayFromBd() const{
    return m_rowsForCellsInner;
  }

  const indices_v_t & graphRowsOfCellsNearBd() const{
    return m_rowsForCellsBd;
  }

private:
  auto boundsImpl(int i) const
  {
    const auto & a = (i==0) ? m_x : (i==1) ? m_y : m_z;
    using numlimits = std::numeric_limits<scalar_type>;
    auto res = std::make_tuple(numlimits::max(), numlimits::min());

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
    pressiodemoapps::impl::readMeshInfo(meshDir, m_dim,
					m_cellDeltas,
					m_cellDeltasInv,
					m_stencilSize,
					m_stencilMeshSize,
					m_sampleMeshSize);

    pressiodemoapps::resize(m_x, m_stencilMeshSize);
    pressiodemoapps::resize(m_y, m_stencilMeshSize);
    pressiodemoapps::resize(m_z, m_stencilMeshSize);
    pressiodemoapps::impl::readMeshCoordinates(meshDir, m_dim, m_x, m_y, m_z);

    const auto graphNumCols = (m_stencilSize-1)*m_dim;
    pressiodemoapps::resize(m_graph, m_sampleMeshSize, graphNumCols+1);
    pressiodemoapps::impl::readMeshConnectivity(meshDir, m_graph, graphNumCols+1);

    // figure out how many cells are near the boundaries
    for (index_t it=0; it<m_sampleMeshSize; ++it)
    {
      if (m_dim==1)
      {
	const auto b1 = hasBdLeft1d(it);
	const auto b2 = hasBdRight1d(it);
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
	const auto b2 = hasBdTop(it);
	const auto b3 = hasBdRight2d(it);
	const auto b4 = hasBdBottom(it);
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

	if(m_stencilSize==7){
	  throw std::runtime_error("3D with 7pt stencil not supported yet.");
	}

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

  bool hasBdLeft1d(const index_t rowInd) const{
    if (m_graph(rowInd, 1)==-1) {
      return true;
    }

    if (m_stencilSize==7){
      if (m_graph(rowInd, 3)==-1) return true;
      if (m_graph(rowInd, 5)==-1) return true;
    }
    return false;
  }

  bool hasBdLeft2d(const index_t rowInd) const{
    if (m_graph(rowInd, 1)==-1) {
      return true;
    }

    if (m_stencilSize==7){
      if (m_graph(rowInd, 5)==-1) return true;
      if (m_graph(rowInd, 9)==-1) return true;
    }
    return false;
  }

  bool hasBdTop(const index_t rowInd) const{
    if (m_graph(rowInd, 2)==-1) {
      return true;
    }

    if (m_stencilSize==7){
      if (m_graph(rowInd, 6)==-1) return true;
      if (m_graph(rowInd, 10)==-1) return true;
    }
    return false;
  }

  bool hasBdRight1d(const index_t rowInd) const{
    if (m_graph(rowInd, 2)==-1) {
      return true;
    }

    if (m_stencilSize==7){
      if (m_graph(rowInd, 4)==-1) return true;
      if (m_graph(rowInd, 6)==-1) return true;
    }
    return false;
  }

  bool hasBdRight2d(const index_t rowInd) const{
    if (m_graph(rowInd, 3)==-1) {
      return true;
    }

    if (m_stencilSize==7){
      if (m_graph(rowInd, 7)==-1) return true;
      if (m_graph(rowInd, 11)==-1) return true;
    }
    return false;
  }

  bool hasBdBottom(const index_t rowInd) const{
    if (m_graph(rowInd, 4)==-1) {
      return true;
    }

    if (m_stencilSize==7){
      if (m_graph(rowInd, 8)==-1) return true;
      if (m_graph(rowInd, 12)==-1) return true;
    }
    return false;
  }

  bool hasBdLeft3d(const index_t rowInd) const{
    if (m_graph(rowInd, 1)==-1) return true;
    return false;
  }

  bool hasBdFront3d(const index_t rowInd) const{
    if (m_graph(rowInd, 2)==-1) return true;
    return false;
  }

  bool hasBdRight3d(const index_t rowInd) const{
    if (m_graph(rowInd, 3)==-1) return true;
    return false;
  }

  bool hasBdBack3d(const index_t rowInd) const{
    if (m_graph(rowInd, 4)==-1) return true;
    return false;
  }

  bool hasBdBottom3d(const index_t rowInd) const{
    if (m_graph(rowInd, 5)==-1) return true;
    return false;
  }

  bool hasBdTop3d(const index_t rowInd) const{
    if (m_graph(rowInd, 6)==-1) return true;
    return false;
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
    graph: contains a list such that

    first col   : contains GIDs of cells where we want velocity
    col 1,2,3,4 : contains GIDs of neighboring cells needed for stencil
    the order of the neighbors is: west, north, east, south
    if needed:
       col 4,5,6,7 : GIDs of degree2 neighboring cells needed for stencil
       the order of the neighbors is: west, north, east, south

       col 8,9,10,11: GIDs of degree3 neighboring cells needed for stencil
       the order of the neighbors is: west, north, east, south
  */
  graph_t m_graph = {};

  indices_v_t m_rowsForCellsInner;
  indices_v_t m_rowsForCellsBd;

};

}}
#endif
