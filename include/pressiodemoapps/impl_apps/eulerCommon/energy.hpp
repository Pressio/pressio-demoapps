
#ifndef PRESSIODEMOAPPS_EE_ENERGY_HPP_
#define PRESSIODEMOAPPS_EE_ENERGY_HPP_

namespace pressiodemoapps{ namespace ee{

/* 1d euler, 3 conserved variables */
template<class scalar_type>
scalar_type computeEnergyFromPrimitive2(const scalar_type & gammaMinusOneInv,
				       const std::array<scalar_type, 3> & prim)
{
  const auto one = static_cast<scalar_type>(1);
  const auto two = static_cast<scalar_type>(2);
  const auto half = one/two;

  const auto usq = prim[1]*prim[1];
  return prim[2]*gammaMinusOneInv + half*prim[0]*(usq);
}

template<class scalar_type>
scalar_type computeEnergyFromPrimitive(const scalar_type & gamma,
				       const std::array<scalar_type, 3> & prim)
{
  const auto one = static_cast<scalar_type>(1);

  const auto gammaMinusOneInv = one/(gamma - one);
  return computeEnergyFromPrimitive2(gammaMinusOneInv, prim);
}

/* 2d euler, 4 conserved variables */
template<class scalar_type>
scalar_type computeEnergyFromPrimitive2(const scalar_type & gammaMinusOneInv,
				       const std::array<scalar_type, 4> & prim)
{
  const auto one = static_cast<scalar_type>(1);
  const auto two = static_cast<scalar_type>(2);
  const auto half = one/two;

  const auto usq = prim[1]*prim[1];
  const auto vsq = prim[2]*prim[2];
  return prim[3]*gammaMinusOneInv + half*prim[0]*(usq+vsq);
}

template<class scalar_type>
scalar_type computeEnergyFromPrimitive(const scalar_type & gamma,
				       const std::array<scalar_type, 4> & prim)
{
  const auto one = static_cast<scalar_type>(1);

  const auto gammaMinusOneInv = one/(gamma - one);
  return computeEnergyFromPrimitive2(gammaMinusOneInv, prim);
}

/* 3d euler, 5 conserved variables */
template<class scalar_type>
scalar_type computeEnergyFromPrimitive2(const scalar_type & gammaMinusOneInv,
				       const std::array<scalar_type, 5> & prim)
{
  const auto one = static_cast<scalar_type>(1);
  const auto two = static_cast<scalar_type>(2);
  const auto half = one/two;

  const auto usq = prim[1]*prim[1];
  const auto vsq = prim[2]*prim[2];
  const auto wsq = prim[3]*prim[3];
  return prim[4]*gammaMinusOneInv + half*prim[0]*(usq+vsq+wsq);
}

template<class scalar_type>
scalar_type computeEnergyFromPrimitive(const scalar_type & gamma,
				       const std::array<scalar_type, 5> & prim)
{
  const auto one = static_cast<scalar_type>(1);

  const auto gammaMinusOneInv = one/(gamma - one);
  return computeEnergyFromPrimitive2(gammaMinusOneInv, prim);
}

}}
#endif
