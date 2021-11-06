
#ifndef PRESSIODAPPS_HAS_SUBSCRIPT_OPERATOR_HPP_
#define PRESSIODAPPS_HAS_SUBSCRIPT_OPERATOR_HPP_

#include <type_traits>

namespace pressiodemoapps{ namespace predicates {

template <typename T, typename arg_t, typename enable = void>
struct is_subscriptable
  : std::false_type{};

template <typename T, typename arg_t>
struct is_subscriptable<
  T, arg_t,
  typename std::enable_if<
    !std::is_void<
      decltype(std::declval< T &>()[std::declval<arg_t const>()])
      >::value and
    !std::is_void<
      decltype(std::declval< T const &>()[std::declval<arg_t const>()])
      >::value
    >::type
  > : std::true_type{};

}}
#endif
