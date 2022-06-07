
#ifndef PRESSIODEMOAPPS_ADVECTION_DIFFUSION_2D_IC_HPP_
#define PRESSIODEMOAPPS_ADVECTION_DIFFUSION_2D_IC_HPP_

namespace pressiodemoapps{ namespace impladvdiff2d{

template<class state_type, class mesh_t, class sc_t>
void burgers2d_gaussian(state_type & state,
			const mesh_t & meshObj,
			sc_t pulseMag,
			sc_t pulseSpread,
			sc_t x0,
			sc_t y0)
{
  // numDofsPerCell = 2 for Burgers2d
  constexpr int numDofPerCell = 2;

  constexpr auto zero = static_cast<sc_t>(0);
  constexpr auto one  = static_cast<sc_t>(1);
  constexpr auto two  = static_cast<sc_t>(2);
  constexpr auto half = one/two;

  const auto &x= meshObj.viewX();
  const auto &y= meshObj.viewY();

  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto ind = i*numDofPerCell;

      const auto dx   = x(i) - x0;
      const auto dy   = y(i) - y0;
      const auto dxSq = dx*dx;
      const auto dySq = dy*dy;
      state(ind)   = pulseMag * std::exp( -(dxSq+dySq)/pulseSpread );
      state(ind+1) = pulseMag * std::exp( -(dxSq+dySq)/pulseSpread );
    }
}

}}//end namespace
#endif
