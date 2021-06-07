
#ifndef PRESSIODEMOAPPS_CCU_MESH_HPP_
#define PRESSIODEMOAPPS_CCU_MESH_HPP_

#include <tuple>
#include "./read_info.hpp"
#include "./read_coords.hpp"
#include "./read_connectivity.hpp"
//#include "./mesh/help_fncs.hpp"

namespace pressiodemoapps{ namespace impl{

template<
  class scalar_type,
  class index_type,
  class xy_type,
  class mesh_g_type,
  bool is_binding
  >
class CellCenteredUniformMesh
{
public:
  using scalar_t  = scalar_type;
  using index_t	  = index_type;
  using x_t	  = xy_type;
  using y_t	  = xy_type;
  using graph_t   = mesh_g_type;

  CellCenteredUniformMesh() = delete;

  template<
    bool _is_binding = is_binding,
    typename std::enable_if<!_is_binding>::type * = nullptr
    >
  explicit CellCenteredUniformMesh(const std::string & meshDir)
  {
    allocateAndSetup(meshDir);
  }

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
  // note that when doing bindings, I need to first construct
  // with {1,1} just so that the numpy array picks up they
  // are 2dim array otherwise it thinks they are 1d array.
  template<
    bool _is_binding = is_binding,
    typename std::enable_if<_is_binding>::type * = nullptr
    >
  explicit CellCenteredUniformMesh(const std::string & meshDir)
    : m_x(1), m_y(1), m_graph({10,5})
  {
    allocateAndSetup(meshDir);
  }
#endif

  void checkStencilSupportsOrder(pressiodemoapps::reconstructionEnum eIn) const
  {
    if (m_stencilSize == 3)
    {
      if (eIn == reconstructionEnum::fifthOrderWeno){
	throw std::runtime_error
	  ("Mesh stencil size not large enough for target reconstruction order.");
      }
    }
  }

  bool isPeriodic() const{
    const auto gSize = (m_stencilSize-1)*m_dim;

    // mesh connectivity is periodi if all neighboring GIDs are positive
    for (index_t i=0; i<m_sampleMeshSize; ++i){
      for (index_t j=0; j<gSize+1; ++j){
	if (m_graph(i, j) < 0){
	  return false;
	}
      }
    }
    return true;
  }

  auto boundsX() const{
    auto res = std::make_tuple(std::numeric_limits<scalar_type>::max(),
			       std::numeric_limits<scalar_type>::min());

    for (int i=0; i<m_stencilMeshSize; ++i){
      auto & v1 = res.get<0>(res);
      auto & v2 = res.get<1>(res);
      v1 = std::min(v1, m_x(i));
      v2 = std::max(v2, m_x(i));
    }
    return res;
  }

  auto boundsY() const{
    auto res = std::make_tuple(std::numeric_limits<scalar_type>::max(),
			       std::numeric_limits<scalar_type>::min());

    for (int i=0; i<m_stencilMeshSize; ++i){
      auto & v1 = res.get<0>(res);
      auto & v2 = res.get<1>(res);
      v1 = std::min(v1, m_y(i));
      v2 = std::max(v2, m_y(i));
    }
    return res;
  }

  int dimensionality() const{
    return m_dim;
  }

  index_t stencilMeshSize() const{
    return m_stencilMeshSize;
  }

  index_t sampleMeshSize() const{
    return m_sampleMeshSize;
  }

  const int stencilSize() const {
    return m_stencilSize;
  }

  const graph_t & graph() const{
    return m_graph;
  }

  const scalar_type dx() const{
    return m_cellDeltas[0];
  }

  const scalar_type dxInv() const{
    return m_cellDeltas[1];
  }

  const scalar_type dy() const{
    return m_cellDeltas[2];
  }

  const scalar_type dyInv() const{
    return m_cellDeltas[3];
  }

  const x_t & viewX() const{
    return m_x;
  }

  const y_t & viewY() const{
    return m_y;
  }

private:
  void allocateAndSetup(const std::string & meshDir)
  {
    pressiodemoapps::impl::readMeshInfo(meshDir, m_dim,
					m_cellDeltas,
					m_stencilSize,
					m_stencilMeshSize,
					m_sampleMeshSize);

    pressiodemoapps::resize(m_x, m_stencilMeshSize);
    pressiodemoapps::resize(m_y, m_stencilMeshSize);
    pressiodemoapps::impl::readMeshCoordinates(meshDir, m_x, m_y);

    const auto graphSize = (m_stencilSize-1)*m_dim;
    pressiodemoapps::resize(m_graph, m_sampleMeshSize, graphSize+1);
    pressiodemoapps::impl::readMeshConnectivity(meshDir, m_graph, graphSize+1);
  }

private:
  int m_dim = {};
  std::array<scalar_type,4> m_cellDeltas{};

  index_t m_stencilMeshSize = {};
  index_t m_sampleMeshSize  = {};

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
  int m_stencilSize = {};

  x_t m_x = {};
  y_t m_y = {};
};

}}
#endif
