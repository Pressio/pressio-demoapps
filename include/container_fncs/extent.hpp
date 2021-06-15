
#ifndef PRESSIODEMOAPPS_EXTENT_FUNC_HPP_
#define PRESSIODEMOAPPS_EXTENT_FUNC_HPP_

namespace pressiodemoapps{

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class value_type>
auto extent(const Eigen::Matrix<value_type, -1, 1> & objIn,
	    int axis)
{
  return (axis==0) ? objIn.size() : 1;
}

template<class value_type>
auto extent(Eigen::Matrix<value_type, -1, 1> & objIn,
	    int axis)
{
  return (axis==0) ? objIn.size() : 1;
}

template<class value_type>
auto extent(const Eigen::Matrix<value_type, -1, -1, Eigen::RowMajor> & objIn,
	    int axis)
{
  return (axis==0) ? objIn.rows() : objIn.cols();
}
template<class value_type>
auto extent(Eigen::Matrix<value_type, -1, -1, Eigen::RowMajor> & objIn,
	    int axis)
{
  return (axis==0) ? objIn.rows() : objIn.cols();
}

template<class value_type>
auto extent(const Eigen::Matrix<value_type, -1, -1, Eigen::ColMajor> & objIn,
	    int axis)
{
  return (axis==0) ? objIn.rows() : objIn.cols();
}
template<class value_type>
auto extent(Eigen::Matrix<value_type, -1, -1, Eigen::ColMajor> & objIn,
	    int axis)
{
  return (axis==0) ? objIn.rows() : objIn.cols();
}
#endif


#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<
  class T,
  typename std::enable_if<
    pressiodemoapps::predicates::is_array_pybind<T>::value
    >::type * = nullptr
  >
auto extent(T & objIn, int axis)
{
  return objIn.shape(axis);
}

template<
  class T,
  typename std::enable_if<
    pressiodemoapps::predicates::is_array_pybind<T>::value
    >::type * = nullptr
  >
auto extent(const T & objIn, int axis)
{
  return objIn.shape(axis);
}

#endif

}
#endif
