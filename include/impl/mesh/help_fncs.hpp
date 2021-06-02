
#ifndef PRESSIODEMOAPPS_MESH_FUNCS_HPP_
#define PRESSIODEMOAPPS_MESH_FUNCS_HPP_

namespace pressiodemoapps{ namespace impl{

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class scalar_type, class sizetype>
void resize(Eigen::Matrix<scalar_type, -1, 1> & objIn,
	    sizetype newSize)
{
  objIn.resize(newSize);
}

template<class int_type, class sizetype>
void resize(Eigen::Matrix<int_type, -1, -1, Eigen::RowMajor> & objIn,
	    sizetype rows,
	    sizetype cols)
{
  objIn.resize(rows, cols);
}
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class T, class sizetype>
typename std::enable_if<
  pressiodemoapps::predicates::is_array_pybind<T>::value
  >::type
resize(T & objIn, sizetype newSize)
{
  if (objIn.ndim() != 1 ){
    throw std::runtime_error
      ("Wrong resize overload called on pybind array");
  }

  objIn.resize({newSize}, false);
}

template<class T, class sizetype>
typename std::enable_if<
  pressiodemoapps::predicates::is_array_pybind<T>::value
  >::type
resize(T & objIn, sizetype rows, sizetype cols)
{
  if (objIn.ndim() != 2 ){
    throw std::runtime_error
      ("Wrong resize overload called on pybind array");
  }
  objIn.resize({rows, cols}, false);
}
#endif


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

  graph_t m_graph;
  int m_loc;

  Neighbors() = delete;

  Neighbors(const mesh_t & meshIn, int loc)
    : m_graph(meshIn.graph()), m_loc(loc)
  {}

  // const auto & operator[](int k) const{
  //   return m_graph(m_loc, k);
  // }

  const auto & operator()(int k) const{
    return m_graph(m_loc, k);
  }
};

template<class mesh_t>
Neighbors<mesh_t> neighbors(const mesh_t & mesh, int loc){
  return Neighbors<mesh_t>(mesh, loc);
}
#endif

}}
#endif
