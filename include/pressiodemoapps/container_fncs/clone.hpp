
#ifndef PRESSIODEMOAPPS_CLONE_FUNC_HPP_
#define PRESSIODEMOAPPS_CLONE_FUNC_HPP_

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace pressiodemoapps{

template<class value_type>
auto clone(const Eigen::Matrix<value_type, Eigen::Dynamic, 1> & o){
  return Eigen::Matrix<value_type, Eigen::Dynamic, 1>(o.size());
}

template<class value_type>
auto clone(const Eigen::Matrix<
	   value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor
	   > & o)
{
  using ret_t = Eigen::Matrix<
    value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
  return ret_t(o.rows(), o.cols());
}

template<class value_type>
auto clone(const Eigen::Matrix<
	   value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & o)
{
  using ret_t = Eigen::Matrix<
    value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  return ret_t(o.rows(), o.cols());
}

template<class value_type, class IntType>
auto clone(const Eigen::SparseMatrix<value_type, Eigen::RowMajor, IntType> & M){
  Eigen::SparseMatrix<value_type, Eigen::RowMajor, IntType> r(M);
  return r;
}

template<class value_type, class IntType>
auto clone(const Eigen::SparseMatrix<value_type, Eigen::ColMajor, IntType> & M){
  Eigen::SparseMatrix<value_type, Eigen::ColMajor, IntType> r(M);
  return r;
}

}
#endif
