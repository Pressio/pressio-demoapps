
#ifndef PRESSIODEMOAPPS_WENO_HPP_
#define PRESSIODEMOAPPS_WENO_HPP_

#include <cmath>
#include <limits>
#include "array"
#include "./impl/weno/weno5.hpp"

namespace pressiodemoapps{

template<typename sc_t>
typename std::enable_if<std::is_floating_point<sc_t>::value>::type
weno5(sc_t & uLeftNeg,
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
  std::is_floating_point<sc_t>::value and
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


template<typename edge_t, typename state_t, typename index_t>
typename std::enable_if<
  !std::is_floating_point<edge_t>::value and
  !std::is_void<decltype(std::declval<state_t>()(0))>::value
  >::type
weno5(edge_t & uLeftNeg,
      edge_t & uLeftPos,
      edge_t & uRightNeg,
      edge_t & uRightPos,
      const state_t & q,
      index_t l2i, index_t l1i, index_t l0i,
      index_t r0i, index_t r1i, index_t r2i,
      index_t uIndex,
      int numDof)
{

  for (int dof=0; dof<numDof; ++dof)
    {
      pressiodemoapps::impl::weno5
	(uLeftNeg(dof), uLeftPos(dof),
	 q(l2i+dof), q(l1i+dof), q(l0i+dof), q(uIndex+dof), q(r0i+dof), q(r1i+dof));

      pressiodemoapps::impl::weno5
	(uRightNeg(dof), uRightPos(dof),
	 q(l1i+dof), q(l0i+dof), q(uIndex+dof), q(r0i+dof), q(r1i+dof), q(r2i+dof));
    }
}

template<typename sc_t, typename state_t, typename index_t>
typename std::enable_if<
  !std::is_void<decltype(std::declval<state_t>()[0])>::value
  >::type
weno5(std::vector<sc_t> & uLeftNeg,
      std::vector<sc_t> & uLeftPos,
      std::vector<sc_t> & uRightNeg,
      std::vector<sc_t> & uRightPos,
      const state_t & q,
      index_t l2i, index_t l1i, index_t l0i,
      index_t r0i, index_t r1i, index_t r2i,
      index_t uIndex,
      int numDof)
{

  for (int dof=0; dof<numDof; ++dof)
    {
      pressiodemoapps::impl::weno5
	(uLeftNeg[dof], uLeftPos[dof],
	 q[l2i+dof], q[l1i+dof], q[l0i+dof], q[uIndex+dof], q[r0i+dof], q[r1i+dof]);

      pressiodemoapps::impl::weno5
	(uRightNeg[dof], uRightPos[dof],
	 q[l1i+dof], q[l0i+dof], q[uIndex+dof], q[r0i+dof], q[r1i+dof], q[r2i+dof]);
    }
}

}
#endif
