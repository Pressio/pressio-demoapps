
#ifndef PRESSIODEMOAPPS_EE_RANKINE_HUGONIOT_HPP_
#define PRESSIODEMOAPPS_EE_RANKINE_HUGONIOT_HPP_

namespace pressiodemoapps{ namespace ee{

// Compute post shock density for a Normal shock wave for a given shock Mach number
template<typename scalar_type>
scalar_type computePostShockDensity(const scalar_type rhoPreShock,
                                    const scalar_type machShock,
                                    const scalar_type gamma)
{
  constexpr auto one = static_cast<scalar_type>(1);
  constexpr auto two = static_cast<scalar_type>(2);
  const auto machShockSq  = machShock * machShock;

  return rhoPreShock * (gamma+one) * machShockSq / (two + (gamma-one) * machShockSq);
}

// Compute post shock pressure for a Normal shock wave for a given shock Mach number
template<typename scalar_type>
scalar_type computePostShockPressure(const scalar_type pPreShock,
                                     const scalar_type machShock,
                                     const scalar_type gamma)
{
  constexpr auto one = static_cast<scalar_type>(1);
  constexpr auto two = static_cast<scalar_type>(2);
  const auto machShockSq  = machShock * machShock;
  return pPreShock * (one + two*gamma/(gamma+one) * (machShockSq-one));
}


// 1D case
template<typename scalar_type>
void computePostShockConditions(std::array<scalar_type,3> & primPostShock,
				const std::array<scalar_type,3> & primPreShock,
				const scalar_type machShock,
				const scalar_type gamma)
{
  const auto rhoPreShock = primPreShock[0];
  const auto velPreShock = primPreShock[1];
  const auto pPreShock   = primPreShock[2];

  const auto rhoPostShock = computePostShockDensity(rhoPreShock,machShock,gamma);
  const auto pPostShock = computePostShockPressure(pPreShock,machShock,gamma);

  // Velocity
  const auto aPreShock = std::sqrt( gamma * pPreShock / rhoPreShock );
  const auto aPostShock = std::sqrt( gamma * pPostShock / rhoPostShock );
  const auto machShockSq  = machShock * machShock;

  constexpr auto one = static_cast<scalar_type>(1);
  constexpr auto two = static_cast<scalar_type>(2);
  constexpr auto half = one/two;

  const auto num = (one + half*(gamma-one)*machShockSq);
  const auto den = gamma*machShockSq - half*(gamma-one);
  const auto machRelPostShock = std::sqrt(num/den);

  // Velocity assuming leftward propogating shock:
  const auto velPostShock = machRelPostShock*aPostShock - machShock*aPreShock + velPreShock;

  // Output to array
  primPostShock[0] = rhoPostShock;
  primPostShock[1] = velPostShock;
  primPostShock[2] = pPostShock;
}


// 2D case
// angle is 0 for a shock propagting along the x axis,
// pi/2 for a shock propagating along the y axis
// angle is wrt to normal from the shock pointing towards the pre-shock region
template<typename scalar_type>
void computePostShockConditionsFromPreshockAtRest(std::array<scalar_type,4> & primPostShock,
						  const std::array<scalar_type,4> & primPreShock,
						  const scalar_type angle,
						  const scalar_type machShock,
						  const scalar_type gamma)
{
  constexpr scalar_type zero{0};
  if (primPreShock[1] != zero and primPreShock[2] != zero){
    throw std::runtime_error
      ("Ranking-Hugoniot conditions functions called with non-at-rest pre shock");
  }

  const auto sinA = std::sin(angle);
  const auto cosA = std::cos(angle);

  std::array<scalar_type, 3> primPreShock1D  = {0., 0., 0.};
  std::array<scalar_type, 3> primPostShock1D = {0., 0., 0.};
  primPreShock1D[0] = primPreShock[0];
  primPreShock1D[1] = zero;
  primPreShock1D[2] = primPreShock[3];
  computePostShockConditions(primPostShock1D, primPreShock1D, machShock, gamma);

  primPostShock[0] = primPostShock1D[0];
  primPostShock[1] = -primPostShock1D[1] * cosA;
  primPostShock[2] = -primPostShock1D[1] * sinA;
  primPostShock[3] = primPostShock1D[2];
}

}}
#endif
