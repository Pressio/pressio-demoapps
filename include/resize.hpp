
#ifndef PRESSIODEMOAPPS_RESIZE_FUNC_HPP_
#define PRESSIODEMOAPPS_RESIZE_FUNC_HPP_

namespace pressiodemoapps{

template<class value_type, class sizetype>
void resize(std::vector<value_type> & objIn,
	    sizetype newSize)
{
  objIn.resize(newSize);
}

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class value_type, class sizetype>
void resize(Eigen::Matrix<value_type, -1, 1> & objIn,
	    sizetype newSize)
{
  objIn.resize(newSize);
}

template<class value_type, class sizetype>
void resize(Eigen::Matrix<value_type, -1, -1, Eigen::RowMajor> & objIn,
	    sizetype rows,
	    sizetype cols)
{
  objIn.resize(rows, cols);
}


template<class value_type, class sizetype>
void resize(Eigen::Matrix<value_type, -1, -1, Eigen::ColMajor> & objIn,
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
      ("Resize overload for ndim=1 called on pybind array with ndim!=1");
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
      ("Resize overload for ndim=2 called on pybind array with ndim!=2");
  }
  objIn.resize({rows, cols}, false);
}
#endif

}
#endif
