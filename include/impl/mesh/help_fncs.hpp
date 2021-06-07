
#ifndef PRESSIODEMOAPPS_MESH_HELP_FUNCS_HPP_
#define PRESSIODEMOAPPS_MESH_HELP_FUNCS_HPP_

namespace pressiodemoapps{ namespace impl{

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class mesh_t>
auto neighbors(const mesh_t & mesh, int loc){
  return mesh.graph().row(loc);
}
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class mesh_t>
struct Neighbors
{
  using index_t = typename mesh_t::index_t;
  using graph_t  = typename mesh_t::mesh_graph_t;

  Neighbors() = delete;
  Neighbors(const mesh_t & meshIn, int loc)
    : m_graph(meshIn.graph()), m_loc(loc)
  {}

  const index_t & operator()(index_t k) const{
    return m_graph(m_loc, k);
  }

private:
  const graph_t & m_graph;
  int m_loc;
};

template<class mesh_t, class T>
Neighbors<mesh_t> neighbors(const mesh_t & mesh, T loc){
  return Neighbors<mesh_t>(mesh, loc);
}
#endif

}}
#endif
