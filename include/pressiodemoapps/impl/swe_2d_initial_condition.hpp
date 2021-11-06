
#ifndef PRESSIODEMOAPPS_SWE2D_IC_HPP_
#define PRESSIODEMOAPPS_SWE2D_IC_HPP_

namespace pressiodemoapps{ namespace implswe{

template<class state_type, class mesh_t, class scalar_type>
void GaussianPulse(state_type & state,
		  const mesh_t & meshObj,
		  const scalar_type pulseMag)
{
  constexpr int numDofPerCell = 3;
  constexpr auto one  = static_cast<scalar_type>(1);

  const auto &x= meshObj.viewX();
  const auto &y= meshObj.viewY();

  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto ind = i*numDofPerCell;
      auto r = std::sqrt( (x(i) - one)*(x(i) - one) + (y(i) - one)*(y(i) - one) );
      state(ind)   = 1. + pulseMag*std::exp( -(r*r) );
      state(ind+1) = 0.;
      state(ind+2) = 0.;
    }
}

}}//end namespace
#endif
