
#ifndef PRESSIODEMOAPPS_NATIVE_PYBIND_ARRAY_HPP_
#define PRESSIODEMOAPPS_NATIVE_PYBIND_ARRAY_HPP_

namespace pressiodemoapps{ namespace predicates {

#include <type_traits>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>

template <typename T, typename enable = void>
struct is_cstyle_array_pybind11
  : std::false_type {};

template <typename T>
struct is_cstyle_array_pybind11<
  T,
  typename std::enable_if<
    std::is_same<
      T,
      pybind11::array_t<
	typename T::value_type,
	pybind11::array::c_style
	>
      >::value
    >::type
  > : std::true_type{};
//----------------------------------------------

template <typename T, typename enable = void>
struct is_fstyle_array_pybind11
  : std::false_type {};

template <typename T>
struct is_fstyle_array_pybind11<
  T,
  typename std::enable_if<
    std::is_same<
      T,
      pybind11::array_t<
	typename T::value_type,
	pybind11::array::f_style
	>
      >::value
    >::type
  > : std::true_type{};
//----------------------------------------------

template <typename T, typename enable = void>
struct is_array_pybind11
  : std::false_type {};

template <typename T>
struct is_array_pybind11<
  T,
  typename std::enable_if<
    is_cstyle_array_pybind11<T>::value or
    is_fstyle_array_pybind11<T>::value
    >::type
  > : std::true_type{};
//----------------------------------------------

template <typename T>
using is_array_pybind = is_array_pybind11<T>;

template <typename T>
using is_cstyle_array_pybind = is_cstyle_array_pybind11<T>;

template <typename T>
using is_fstyle_array_pybind = is_fstyle_array_pybind11<T>;

}}
#endif
