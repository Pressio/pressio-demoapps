
#ifndef PRESSIODEMOAPPS_CCU_MESH_HPP_
#define PRESSIODEMOAPPS_CCU_MESH_HPP_

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
  using scalar_t     = scalar_type;
  using index_t	     = index_type;
  using x_t	     = xy_type;
  using y_t	     = xy_type;
  using mesh_graph_t = mesh_g_type;

  CellCenteredUniformMesh() = delete;

  template<
    bool _is_binding = is_binding,
    typename std::enable_if<!_is_binding>::type * = nullptr
    >
  explicit CellCenteredUniformMesh(const std::string & meshDir)
  {
    resizeAndSetup(meshDir);
  }

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
  template<
    bool _is_binding = is_binding,
    typename std::enable_if<_is_binding>::type * = nullptr
    >
  explicit CellCenteredUniformMesh(const std::string & meshDir)
    : m_x(1), m_y(1), m_graph({1,1})
  {
    resizeAndSetup(meshDir);
  }
#endif

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

  const mesh_graph_t & graph() const{
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

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN

//   // ptIndex is not the global id
//   auto neighborsOf(const index_t smPt) const{
//     return m_graph.row(smPt);
//   }

// #elif defined PRESSIODEMOAPPS_ENABLE_BINDINGS

//   auto neighborsOf(const index_t smPt) const
//   {
//     auto res = m_graph[pybind11::make_tuple(smPt, pybind11::slice(0, m_graph.shape(1), 1))];
//     // auto tp = pybind11::make_tuple(pybind11::slice(0, m_graph.shape(1), 1));
//     // pybind11::array Aslice = m_graph(smPt, tp);
//     return res;
//   }
// #endif

private:
  void resizeAndSetup(const std::string & meshDir)
  {
    pressiodemoapps::impl::readMeshInfo(meshDir, m_dim,
					m_cellDeltas,
					m_stencilSize,
					m_stencilMeshSize,
					m_sampleMeshSize);

    pressiodemoapps::impl::resize(m_x, m_stencilMeshSize);
    pressiodemoapps::impl::resize(m_y, m_stencilMeshSize);
    pressiodemoapps::impl::readMeshCoordinates(meshDir, m_x, m_y);

    const auto graphSize = (m_stencilSize-1)*m_dim;
    pressiodemoapps::impl::resize(m_graph, m_sampleMeshSize, graphSize+1);
    pressiodemoapps::impl::readMeshConnectivity(meshDir, m_graph, graphSize+1);
  }

private:
  int m_dim = {};
  std::array<scalar_type,4> m_cellDeltas{};

  index_t m_stencilMeshSize = {};
  index_t m_sampleMeshSize  = {};

  /*
    graph: contains a list such that
    1 0 3 2 -1

    first col   : contains GIDs of cells where we want velocity
    col 1,2,3,4 : contains GIDs of neighboring cells needed for stencil
    the order of the neighbors is: west, north, east, south
    if needed:
       col 4,5,6,7 : GIDs of degree2 neighboring cells needed for stencil
       the order of the neighbors is: west, north, east, south

       col 8,9,10,11: GIDs of degree3 neighboring cells needed for stencil
       the order of the neighbors is: west, north, east, south
  */
  mesh_graph_t m_graph = {};
  int m_stencilSize = {};

  x_t m_x = {};
  y_t m_y = {};
};

}}
#endif
