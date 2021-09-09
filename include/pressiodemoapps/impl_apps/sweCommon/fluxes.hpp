
#ifndef PRESSIODEMOAPPS_SWE_FLUXES_HPP_
#define PRESSIODEMOAPPS_SWE_FLUXES_HPP_

namespace pressiodemoapps{ namespace swe{ namespace impl{

template<typename T, class normal_t, class sc_t>
void sweRusanovFluxThreeDof(T & F,
			  const T & qL,
			  const T & qR,
			  const normal_t & n,
			  const sc_t gravity)
{
  constexpr auto one  = static_cast<sc_t>(1);
  constexpr auto two  = static_cast<sc_t>(2);
  constexpr auto half = one/two;

  constexpr sc_t es{1.e-30};
  sc_t FL[3];
  sc_t FR[3];

  // left
  const auto hL = qL(0);
  const auto uL = qL(1)/(hL + es);
  const auto vL = qL(2)/(hL + es);
  const auto unL = uL*n[0] + vL*n[1];
  const auto pL = 0.5*gravity*hL*hL;
  FL[0] = hL*unL;
  FL[1] = hL*unL*uL + pL*n[0];
  FL[2] = hL*unL*vL + pL*n[1];

  // right
  const auto hR = qR(0);
  const auto uR = qR(1)/(hR + es);
  const auto vR = qR(2)/(hR + es);
  const auto unR = uR*n[0] + vR*n[1];
  const auto pR = 0.5*gravity*hR*hR;
  FR[0] = hR*unR;
  FR[1] = hR*unR*uR + pR*n[0];
  FR[2] = hR*unR*vR + pR*n[1];

  // compute flux
  const auto hm = 0.5*(hL + hR);
  const auto um = (unL*std::sqrt(hL) + unR*std::sqrt(hR) )/(std::sqrt(hL) + std::sqrt(hR)  + es);
  const auto smax = std::abs(um) + std::abs( std::sqrt(gravity*hm) );

  F(0) = half*( FL[0] + FR[0] + smax*(qL(0) - qR(0)) );
  F(1) = half*( FL[1] + FR[1] + smax*(qL(1) - qR(1)) );
  F(2) = half*( FL[2] + FR[2] + smax*(qL(2) - qR(2)) );
}


}}}
#endif
