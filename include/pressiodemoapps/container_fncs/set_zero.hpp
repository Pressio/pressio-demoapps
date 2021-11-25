
#ifndef PRESSIODEMOAPPS_SET_ZERO_HPP_
#define PRESSIODEMOAPPS_SET_ZERO_HPP_

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace pressiodemoapps{

template<class value_type, int I0, int I1>
void set_zero(Eigen::Matrix<value_type, I0, I1, Eigen::ColMajor> & a){
  a.setZero();
}

template<class value_type, int I0, int I1>
void set_zero(Eigen::Matrix<value_type, I0, I1, Eigen::RowMajor> & a){
  a.setZero();
}

template<class value_type, class IntType>
void set_zero(Eigen::SparseMatrix<value_type, Eigen::RowMajor, IntType> & M)
{
  auto values = M.valuePtr();
  for (int i=0; i<M.nonZeros(); ++i){
    values[i] = static_cast<value_type>(0);
  }
}

template<class value_type, class IntType>
void set_zero(Eigen::SparseMatrix<value_type, Eigen::ColMajor, IntType> & M)
{
  auto values = M.valuePtr();
  for (int i=0; i<M.nonZeros(); ++i){
    values[i] = static_cast<value_type>(0);
  }
}

template<class value_type>
void set_zero(Eigen::Ref<Eigen::Matrix<value_type, -1, 1>> & a)
{
  for (int i=0; i<a.rows(); ++i){
    a(i) = static_cast<value_type>(0);
  }
}

template<class value_type>
void set_zero(Eigen::Ref<Eigen::Matrix<value_type, -1, -1, Eigen::ColMajor>> & a)
{
  for (int i=0; i<a.rows(); ++i){
    for (int j=0; j<a.cols(); ++j){
      a(i,j) = static_cast<value_type>(0);
    }
  }
}

template<class value_type>
void set_zero(Eigen::Ref<Eigen::Matrix<value_type, -1, -1, Eigen::RowMajor>> & a)
{
  for (int i=0; i<a.rows(); ++i){
    for (int j=0; j<a.cols(); ++j){
      a(i,j) = static_cast<value_type>(0);
    }
  }
}

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
template<class value_type, class IntType>
void set_zero_omp(Eigen::SparseMatrix<value_type, Eigen::RowMajor, IntType> & M)
{
  auto values = M.valuePtr();
#pragma omp for schedule(auto)
  for (int i=0; i<M.nonZeros(); ++i){
    values[i] = static_cast<value_type>(0);
  }
}

template<class value_type, class IntType>
void set_zero_omp(Eigen::SparseMatrix<value_type, Eigen::ColMajor, IntType> & M)
{
  auto values = M.valuePtr();
#pragma omp for schedule(auto)
  for (int i=0; i<M.nonZeros(); ++i){
    values[i] = static_cast<value_type>(0);
  }
}

template<class value_type>
void set_zero_omp(Eigen::Matrix<value_type, -1, 1> & V)
{
#pragma omp for schedule(auto)
   for (decltype(V.size()) i=0; i<V.size(); ++i){
     V(i) = static_cast<value_type>(0);
   }
}

template<class value_type>
void set_zero_omp(Eigen::Ref<Eigen::Matrix<value_type, -1, 1>> & V)
{
#pragma omp for schedule(auto)
   for (decltype(V.size()) i=0; i<V.size(); ++i){
     V(i) = static_cast<value_type>(0);
   }
}
#endif

}
#endif
