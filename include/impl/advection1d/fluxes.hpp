
#ifndef PRESSIODEMOAPPS_EE_FLUXES_HPP_
#define PRESSIODEMOAPPS_EE_FLUXES_HPP_

namespace pressiodemoapps{ namespace ad{ namespace impl{

template<typename sc_t>
void linAdvRusanovFlux(sc_t & F,
		       const sc_t & qL,
		       const sc_t & qR)
{
  constexpr sc_t es = 1.e-30;
  sc_t FL;
  sc_t FR;

  FL = qL;
  FR = qR;
  F = 0.5*( FL + FR + (qL - qR) );
  //F = qL;
}

}}}
#endif
