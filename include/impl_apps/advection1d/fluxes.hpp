
#ifndef PRESSIODEMOAPPS_LINDADV_FLUXES_HPP_
#define PRESSIODEMOAPPS_LINDADV_FLUXES_HPP_

namespace pressiodemoapps{ namespace ad{ namespace impl{

template<typename sc_t>
void linAdvRusanovFlux(sc_t & F,
		       const sc_t & qL,
		       const sc_t & qR)
{
  //constexpr auto oneHalf = static_cast<sc_t>(0.5);
  // sc_t FL;
  // sc_t FR;
  // FL = qL;
  // FR = qR;
  F = qL; // = 0.5*( qL + qR + (qL - qR) );
}

}}}
#endif
