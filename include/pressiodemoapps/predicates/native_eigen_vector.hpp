
#ifndef PRESSIODEMOAPPS_PREDICATES_NATIVE_EIGEN_VECTOR_HPP_
#define PRESSIODEMOAPPS_PREDICATES_NATIVE_EIGEN_VECTOR_HPP_

#include <Eigen/Dense>

namespace pressiodemoapps{

template <typename T, typename enable = void>
struct is_dynamic_row_vector_eigen : std::false_type {};

template <typename T>
struct is_dynamic_row_vector_eigen<
  T,
  std::enable_if_t<
    std::is_same<
      typename std::remove_cv<T>::type,
		  Eigen::Matrix<typename T::Scalar,1, Eigen::Dynamic>
		 >::value
    >
  > : std::true_type{};
//----------------------------------------------

template <typename T, typename enable = void>
struct is_static_row_vector_eigen : std::false_type {};

template <typename T>
struct is_static_row_vector_eigen<
  T,
  std::enable_if_t<
    std::is_same<
     typename std::remove_cv<T>::type,
		 Eigen::Matrix<typename T::Scalar, 1, T::ColsAtCompileTime>
		 >::value and
    !is_dynamic_row_vector_eigen<T>::value
    >
  > : std::true_type{};
//----------------------------------------------

template <typename T, typename enable = void>
struct is_dynamic_column_vector_eigen : std::false_type {};

template <typename T>
struct is_dynamic_column_vector_eigen<
  T,
  std::enable_if_t<
    std::is_same<
     typename std::remove_cv<T>::type,
		 Eigen::Matrix<typename T::Scalar, Eigen::Dynamic, 1>
		 >::value
    >
  > : std::true_type{};
//----------------------------------------------

template <typename T, typename enable = void>
struct is_static_column_vector_eigen : std::false_type {};

template <typename T>
struct is_static_column_vector_eigen<
  T,
  std::enable_if_t<
    std::is_same<
     typename std::remove_cv<T>::type,
		 Eigen::Matrix<typename T::Scalar, T::RowsAtCompileTime,1>
		 >::value and
    !is_dynamic_column_vector_eigen<T>::value
    >
  > : std::true_type{};
//----------------------------------------------

template <typename T, typename enable = void>
struct is_static_vector_eigen : std::false_type {};

template <typename T>
struct is_static_vector_eigen<
  T,
  std::enable_if_t<
    is_static_row_vector_eigen<T>::value ||
    is_static_column_vector_eigen<T>::value
    >
  > : std::true_type{};
//----------------------------------------------

template <typename T, typename enable = void>
struct is_dynamic_vector_eigen : std::false_type {};

template <typename T>
struct is_dynamic_vector_eigen<
  T,
  std::enable_if_t<
    is_dynamic_row_vector_eigen<T>::value or
    is_dynamic_column_vector_eigen<T>::value
    >
  > : std::true_type{};
//----------------------------------------------

template <typename T, typename enable = void>
struct is_vector_eigen : std::false_type {};

template <typename T>
struct is_vector_eigen<
  T,
  std::enable_if_t<
    is_dynamic_vector_eigen<T>::value or
    is_static_vector_eigen<T>::value
    >
  > : std::true_type{};

}//end namespace
#endif
