
#ifndef PRESSIODEMOAPPS_RESIZE_FUNC_HPP_
#define PRESSIODEMOAPPS_RESIZE_FUNC_HPP_

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace pressiodemoapps{

template<class value_type, class sizetype>
void resize(std::vector<value_type> & a,
	    sizetype newSize)
{
  a.resize(newSize);
}

template<class value_type, class sizetype>
void resize(Eigen::Matrix<value_type, -1, 1> & a,
	    sizetype newSize)
{
  a.resize(newSize);
}

template<class value_type, class sizetype>
void resize(Eigen::Matrix<value_type, -1, -1, Eigen::RowMajor> & M,
	    sizetype rows,
	    sizetype cols)
{
  M.resize(rows, cols);
}

template<class value_type, class sizetype>
void resize(Eigen::Matrix<value_type, -1, -1, Eigen::ColMajor> & M,
	    sizetype rows,
	    sizetype cols)
{
  M.resize(rows, cols);
}

template<class value_type, class IntType, class sizetype>
void resize(Eigen::SparseMatrix<value_type, Eigen::RowMajor, IntType> & M,
	    sizetype rows,
	    sizetype cols)
{
  M.resize(rows, cols);
}

template<class value_type, class IntType, class sizetype>
void resize(Eigen::SparseMatrix<value_type, Eigen::ColMajor, IntType> & M,
	    sizetype rows,
	    sizetype cols)
{
  M.resize(rows, cols);
}

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class T, class sizetype>
typename std::enable_if<
  pressiodemoapps::predicates::is_array_pybind<T>::value
  >::type
resize(T & o, sizetype newSize)
{
  if (o.ndim() != 1 ){
    throw std::runtime_error
      ("Resize overload for ndim=1 called on pybind array with ndim!=1");
  }

  o.resize({newSize}, false);
}

template<class T, class sizetype>
typename std::enable_if<
  pressiodemoapps::predicates::is_array_pybind<T>::value
  >::type
resize(T & o, sizetype rows, sizetype cols)
{
  if (o.ndim() != 2 ){
    throw std::runtime_error
      ("Resize overload for ndim=2 called on pybind array with ndim!=2");
  }
  o.resize({rows, cols}, false);
}
#endif

}
#endif
