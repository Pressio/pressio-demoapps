
#ifndef PRESSIODEOAPPS_PREDICATES_EIGEN_SPARSE_MATRIX_HPP_
#define PRESSIODEOAPPS_PREDICATES_EIGEN_SPARSE_MATRIX_HPP_

#include "Eigen/Sparse"

namespace pressiodemoapps{

template <typename T, typename enable = void>
struct is_sparse_matrix_eigen : std::false_type {};

/*
 * T is an eigen sparse matrix if is
 * not an eigen vector
 * is same type as sparse matrix
*/
template<typename T>
struct is_sparse_matrix_eigen<
  T,
  std::enable_if_t<
    !is_vector_eigen<T>::value and
    std::is_same<
      typename std::remove_cv<T>::type,
      Eigen::SparseMatrix<typename T::Scalar, T::Options, typename T::StorageIndex>
      >::value
    >
  > : std::true_type{};


//----------------------------------------------------------------------

template <typename T1, typename T2, typename enable = void>
struct sparse_sharedmem_eigen_same_storage : std::false_type{};

template <typename T1, typename T2>
struct sparse_sharedmem_eigen_same_storage<
  T1, T2,
  std::enable_if_t<
    (T1::is_row_major && T2::is_row_major) ||
    (T1::is_col_major && T2::is_col_major)
    >
  > : std::true_type{};

}//end namespace
#endif  // TYPE_TRAITS_NATIVE_EIGEN_SPARSE_MATRIX_HPP_
