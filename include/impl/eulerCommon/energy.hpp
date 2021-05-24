
#ifndef PRESSIODEMOAPPS_EE_ENERGY_HPP_
#define PRESSIODEMOAPPS_EE_ENERGY_HPP_

namespace pressiodemoapps{ namespace ee{

template<class scalar_type>
scalar_type computeEnergyFromPrimitive2(const scalar_type & gammaMinusOneInv,
				       const std::array<scalar_type, 4> & prim)
{
  const auto usq = prim[1]*prim[1];
  const auto vsq = prim[2]*prim[2];
  return prim[3]*gammaMinusOneInv + 0.5*prim[0]*(usq+vsq);
}

template<class scalar_type>
scalar_type computeEnergyFromPrimitive(const scalar_type & gamma,
				       const std::array<scalar_type, 4> & prim)
{
  const scalar_type gammaMinusOneInv = 1./(gamma -1.);
  return computeEnergyFromPrimitive2(gammaMinusOneInv, prim);
}

template<class scalar_type>
scalar_type computeEnergyFromPrimitive2(const scalar_type & gammaMinusOneInv,
				       const std::array<scalar_type, 3> & prim)
{
  const auto usq = prim[1]*prim[1];
  return prim[2]*gammaMinusOneInv + 0.5*prim[0]*(usq);
}

template<class scalar_type>
scalar_type computeEnergyFromPrimitive(const scalar_type & gamma,
				       const std::array<scalar_type, 3> & prim)
{
  const scalar_type gammaMinusOneInv = 1./(gamma -1.);
  return computeEnergyFromPrimitive2(gammaMinusOneInv, prim);
}

}}
#endif
