
#ifndef PRESSIODEMOAPPS_ADV_DIFF_2D_RUSANOV_FLUX_VALUES_FUNC_HPP_
#define PRESSIODEMOAPPS_ADV_DIFF_2D_RUSANOV_FLUX_VALUES_FUNC_HPP_

namespace pressiodemoapps{ namespace impladvdiff2d{

template<class sc_t, typename T>
void burgers_rusanov_flux_2d(T & F,
			     const T & qL,
			     const T & qR)
{

  const sc_t one  = static_cast<sc_t>(1);
  const sc_t two  = static_cast<sc_t>(2);
  const sc_t four = static_cast<sc_t>(4);
  const sc_t fourInv = one/four;
  const sc_t lambda = std::abs(qL(0)+qR(0))*fourInv;
  F(0) = fourInv*(qL(0)*qL(0) + qR(0)*qR(0)) + lambda*(qL(0) - qR(0));
}

template<class sc_t, typename T1, class T2>
void burgers_rusanov_flux_jacobian_2d(T1 & JL, T1 & JR,
				      const T2 & qL,
				      const T2 & qR)
{
  const sc_t one  = static_cast<sc_t>(1);
  const sc_t two  = static_cast<sc_t>(2);
  const sc_t four = static_cast<sc_t>(4);
  const sc_t twoInv = one/two;
  const sc_t fourInv = one/four;

  const sc_t qLqRDiff = qL(0) - qR(0);
  const sc_t lambda = std::abs(qL(0)+qR(0));
  const sc_t dLamb = (qL(0)+qR(0)) / lambda;

  JL(0,0) = twoInv * qL(0) + fourInv * dLamb * qLqRDiff + fourInv*lambda;
  JR(0,0) = twoInv * qR(0) + fourInv * dLamb * qLqRDiff - fourInv*lambda;
}

}}
#endif
