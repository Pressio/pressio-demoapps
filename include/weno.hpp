
#ifndef PRESSIODEMOAPPS_WENO_HPP_
#define PRESSIODEMOAPPS_WENO_HPP_

#include <cmath>
#include <limits>
#include "array"

#include "./impl/weno/weno5.hpp"

namespace pressiodemoapps{

template<typename sc_t>
void weno5(sc_t & uLeftNeg,
	   sc_t & uLeftPos,
	   sc_t & uRightNeg,
	   sc_t & uRightPos,
	   const std::array<sc_t, 7> & q)
{
  pressiodemoapps::impl::weno5
    (uLeftNeg, uLeftPos,
     q[0], q[1], q[2], q[3], q[4], q[5]);

  pressiodemoapps::impl::weno5
    (uRightNeg, uRightPos,
     q[1], q[2], q[3], q[4], q[5], q[6]);
}

template<typename sc_t, typename state_t>
typename std::enable_if<
  !std::is_void<decltype(std::declval<state_t>()(0))>::value
  >::type
weno5(sc_t & uLeftNeg,
      sc_t & uLeftPos,
      sc_t & uRightNeg,
      sc_t & uRightPos,
      const state_t & q)
{
  pressiodemoapps::impl::weno5
    (uLeftNeg, uLeftPos,
     q(0), q(1), q(2), q(3), q(4), q(5));

  pressiodemoapps::impl::weno5
    (uRightNeg, uRightPos,
     q(1), q(2), q(3), q(4), q(5), q(6));
}

// template<typename sc_t, typename state_t, typename index_t>
// void weno5(sc_t & uLeftNeg,
//       sc_t & uLeftPos,
//       sc_t & uRightNeg,
//       sc_t & uRightPos,
//       const state_t & q,
//       index_t ind)
// {
//   pressiodemoapps::impl::weno5
//     (uLeftNeg, uLeftPos,
//      q[ind-3], q[ind-2], q[ind-1], q[ind], q[ind+1], q[ind+2]);

//   pressiodemoapps::impl::weno5
//     (uRightNeg, uRightPos,
//      q[ind-2], q[ind-1], q[ind], q[ind+1], q[ind+2], q[ind+3]);
// }

// template<typename sc_t, typename state_t, typename index_t>
// typename std::enable_if<
//   !std::is_void<decltype(std::declval<state_t>()(0))>::value
//   // ::pressiodemoapps::predicates::is_subscriptable<state_t>::value
//   >::type
// weno5(sc_t & uLeftNeg,
//       sc_t & uLeftPos,
//       sc_t & uRightNeg,
//       sc_t & uRightPos,
//       const state_t & q,
//       index_t ind)
// {
//   pressiodemoapps::impl::weno5
//     (uLeftNeg, uLeftPos,
//      q(ind-3), q(ind-2), q(ind-1), q(ind), q(ind+1), q(ind+2));

//   pressiodemoapps::impl::weno5
//     (uRightNeg, uRightPos,
//      q(ind-2), q(ind-1), q(ind), q(ind+1), q(ind+2), q(ind+3));
// }

}
#endif
