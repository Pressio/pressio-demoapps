
#ifndef PRESSIODEMOAPPS_WENO3_IMPL_HPP_
#define PRESSIODEMOAPPS_WENO3_IMPL_HPP_

#include <cmath>

namespace pressiodemoapps{ namespace impl{

template<typename sc_t>
void weno3(sc_t & uNeg,
	   sc_t & uPos,
	   const sc_t & qim1,
	   const sc_t & qi,
	   const sc_t & qip1,
	   const sc_t & qip2)
{

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // ENO stencils:
  //
  //           |______S1_______|
  //           |               |
  //   |_______S0______|       |
  // ..|---o---|---o---|---o---|
  //   | I{i-1}|  I{i} | I{i+1}|
  //                  -|
  //                  i+1/2
  //
  //           |_______S0______|
  //           |               |
  //           |       |_______S1______|
  //         ..|---o---|---o---|...
  //           |  I{i} | I{i+1}| I{i+2}|
  //                   |+
  //                  i+1/2
  //
  // WENO stencil: S{i} = [ I{i-1}, ..., I{i+2} ]
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // see e.g. https://github.com/jbnunezd/fd-solver-mhd/blob/main/src/reconstruction.f90
  // see https://www3.nd.edu/~zxu2/acms60790S13/Shu-WENO-notes.pdf

  constexpr auto epsilon  = static_cast<sc_t>(1e-6);
  constexpr auto one	  = static_cast<sc_t>(1);
  constexpr auto two	  = static_cast<sc_t>(2);
  constexpr auto three	  = static_cast<sc_t>(3);
  constexpr auto oneOvtwo   = one/two;
  constexpr auto oneOvthree = one/three;
  constexpr auto twoOvthree = two/three;

  {
    // u_{i+1/2}^{-}
    const auto p0 = (-qim1 + three*qi)*oneOvtwo;
    const auto p1 = (qi   + qip1)*oneOvtwo;

    const auto B0 = std::pow(qim1-qi, two);
    const auto B1 = std::pow(qi-qip1, two);

    const auto alpha0    = oneOvthree/std::pow(epsilon + B0, two);
    const auto alpha1    = twoOvthree/std::pow(epsilon + B1, two);
    const auto alphaSInv = one/(alpha0 + alpha1);

    const auto w0 = alpha0*alphaSInv;
    const auto w1 = alpha1*alphaSInv;
    uNeg = w0*p0 + w1*p1;
  }

  {
    // u_{i+1/2}^{+}
    const auto p0 = ( qi + qip1)*oneOvtwo;
    const auto p1 = ( three*qip1 - qip2)*oneOvtwo;

    const auto B0 = std::pow(qi-qip1, two);
    const auto B1 = std::pow(qip1-qip2, two);

    const auto alpha0    = twoOvthree/std::pow(epsilon + B0, two);
    const auto alpha1    = oneOvthree/std::pow(epsilon + B1, two);
    const auto alphaSInv = one/(alpha0 + alpha1);

    const auto w0 = alpha0*alphaSInv;
    const auto w1 = alpha1*alphaSInv;
    uPos = w0*p0 + w1*p1;
  }

}

}}
#endif
