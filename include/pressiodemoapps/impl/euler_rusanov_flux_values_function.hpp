
#ifndef PRESSIODEMOAPPS_EE_FLUXES_HPP_
#define PRESSIODEMOAPPS_EE_FLUXES_HPP_

namespace pressiodemoapps{ namespace ee{ namespace impl{

template<class T, typename sc_t>
void eeRusanovFluxThreeDof(T & F,
			   const T & qL,
			   const T & qR,
			   const sc_t gamma)
{
  constexpr auto one  = static_cast<sc_t>(1);
  constexpr auto two  = static_cast<sc_t>(2);
  constexpr auto half = one/two;

  constexpr sc_t es{1.e-30};
  sc_t FL[3];
  sc_t FR[3];

  // left
  const auto rL = qL(0);
  const auto uL = qL(1)/(rL + es);
  const auto pL = (gamma-1)*( qL(2) - half*rL*(uL*uL) );
  const auto HL = ( qL(2) + pL ) / rL;
  FL[0] = rL*uL;
  FL[1] = rL*uL*uL + pL;
  FL[2] = rL*uL*HL;

  // right
  const auto rR = qR(0);
  const auto uR = qR(1)/(rR + es);
  const auto pR = (gamma-1)*( qR(2) - half*rR*(uR*uR) );
  const auto HR = ( qR(2) + pR ) / rR;
  FR[0] = rR*uR;
  FR[1] = rR*uR*uR + pR;
  FR[2] = rR*uR*HR;

  // compute flux
  const auto RT = std::sqrt(rR/(rL));
  const auto u = (uL+RT*uR)/(1.+RT);
  const auto H = ( HL+RT* HR)/(1.+RT);
  const auto a = std::sqrt( (gamma-1.)*(H - half*(u*u)) );
  const auto smax = std::abs(u) + a;

  F(0) = half*( FL[0] + FR[0] + smax*(qL(0) - qR(0)) );
  F(1) = half*( FL[1] + FR[1] + smax*(qL(1) - qR(1)) );
  F(2) = half*( FL[2] + FR[2] + smax*(qL(2) - qR(2)) );
}

template<typename T, class normal_t, class sc_t>
void eeRusanovFluxFourDof(T & F,
			  const T & qL,
			  const T & qR,
			  const normal_t & n,
			  const sc_t gamma)
{
  constexpr auto one  = static_cast<sc_t>(1);
  constexpr auto two  = static_cast<sc_t>(2);
  constexpr auto half = one/two;

  constexpr sc_t es{1.e-30};
  sc_t FL[4];
  sc_t FR[4];

  // left
  const auto rL = qL(0);
  const auto uL = qL(1)/(rL + es);
  const auto vL = qL(2)/(rL + es);
  const auto unL = uL*n[0] + vL*n[1];
  const auto pL = (gamma-1)*( qL(3) - half*rL*(uL*uL+vL*vL) );
  const auto HL = ( qL(3) + pL ) / rL;
  FL[0] = rL*unL;
  FL[1] = rL*unL*uL + pL*n[0];
  FL[2] = rL*unL*vL + pL*n[1];
  FL[3] = rL*unL*HL;

  // right
  const auto rR = qR(0);
  const auto uR = qR(1)/(rR + es);
  const auto vR = qR(2)/(rR + es);
  const auto unR = uR*n[0] + vR*n[1];
  const auto pR = (gamma-1)*( qR(3) - half*rR*(uR*uR+vR*vR) );
  const auto HR = ( qR(3) + pR ) / rR;
  FR[0] = rR*unR;
  FR[1] = rR*unR*uR + pR*n[0];
  FR[2] = rR*unR*vR + pR*n[1];
  FR[3] = rR*unR*HR;

  // compute flux
  const auto RT = std::sqrt(rR/(rL));
  const auto u = (uL+RT*uR)/(1.+RT);
  const auto v = (vL+RT*vR)/(1.+RT);
  const auto H = ( HL+RT* HR)/(1.+RT);
  const auto a = std::sqrt( (gamma-1.)*(H - half*(u*u+v*v)) );
  const auto smax = std::sqrt(u*u+v*v) + a;

  F(0) = half*( FL[0] + FR[0] + smax*(qL(0) - qR(0)) );
  F(1) = half*( FL[1] + FR[1] + smax*(qL(1) - qR(1)) );
  F(2) = half*( FL[2] + FR[2] + smax*(qL(2) - qR(2)) );
  F(3) = half*( FL[3] + FR[3] + smax*(qL(3) - qR(3)) );
}

template<typename T, class normal_t, class sc_t>
void eeRusanovFluxFiveDof(T & F,
			  const T & qL,
			  const T & qR,
			  const normal_t & n,
			  const sc_t gamma)
{
  constexpr auto one  = static_cast<sc_t>(1);
  constexpr auto two  = static_cast<sc_t>(2);
  constexpr auto half = one/two;

  constexpr sc_t es{1.e-30};
  sc_t FL[5];
  sc_t FR[5];

  // left
  const auto rL = qL(0);
  const auto uL = qL(1)/(rL + es);
  const auto vL = qL(2)/(rL + es);
  const auto wL = qL(3)/(rL + es);
  const auto unL = uL*n[0] + vL*n[1] + wL*n[2];
  const auto pL = (gamma-1)*( qL(4) - half*rL*(uL*uL+vL*vL+wL*wL) );
  const auto HL = ( qL(4) + pL ) / rL;
  FL[0] = rL*unL;
  FL[1] = rL*unL*uL + pL*n[0];
  FL[2] = rL*unL*vL + pL*n[1];
  FL[3] = rL*unL*wL + pL*n[2];
  FL[4] = rL*unL*HL;

  // right
  const auto rR = qR(0);
  const auto uR = qR(1)/(rR + es);
  const auto vR = qR(2)/(rR + es);
  const auto wR = qR(3)/(rR + es);
  const auto unR = uR*n[0] + vR*n[1] + wR*n[2];
  const auto pR = (gamma-1)*( qR(4) - half*rR*(uR*uR+vR*vR+wR*wR) );
  const auto HR = ( qR(4) + pR ) / rR;
  FR[0] = rR*unR;
  FR[1] = rR*unR*uR + pR*n[0];
  FR[2] = rR*unR*vR + pR*n[1];
  FR[3] = rR*unR*wR + pR*n[2];
  FR[4] = rR*unR*HR;

  // compute flux
  const auto RT = std::sqrt(rR/(rL));
  const auto u = (uL+RT*uR)/(1.+RT);
  const auto v = (vL+RT*vR)/(1.+RT);
  const auto w = (wL+RT*wR)/(1.+RT);
  const auto H = ( HL+RT* HR)/(1.+RT);
  const auto a = std::sqrt( (gamma-1.)*(H - half*(u*u + v*v + w*w)) );
  const auto smax = std::sqrt(u*u+v*v+w*w) + a;

  F(0) = half*( FL[0] + FR[0] + smax*(qL(0) - qR(0)) );
  F(1) = half*( FL[1] + FR[1] + smax*(qL(1) - qR(1)) );
  F(2) = half*( FL[2] + FR[2] + smax*(qL(2) - qR(2)) );
  F(3) = half*( FL[3] + FR[3] + smax*(qL(3) - qR(3)) );
  F(4) = half*( FL[4] + FR[4] + smax*(qL(4) - qR(4)) );
}

}}}
#endif
