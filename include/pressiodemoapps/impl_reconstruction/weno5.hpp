
#ifndef PRESSIODEMOAPPS_WENO5_IMPL_HPP_
#define PRESSIODEMOAPPS_WENO5_IMPL_HPP_

#include <cmath>

namespace pressiodemoapps{ namespace impl{

template<typename sc_t>
void weno5(sc_t & uNeg,
	   sc_t & uPos,
	   const sc_t & qim2,
	   const sc_t & qim1,
	   const sc_t & qi,
	   const sc_t & qip1,
	   const sc_t & qip2,
	   const sc_t & qip3)
{

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  // 	 |           |   u(i)    |           |
  // 	 |  u(i-1)   |___________|           |
  // 	 |___________|           |   u(i+1)  |
  // 	 |           |           |___________|
  // ...|-----0-----|-----0-----|-----0-----|...
  // 	 |    i-1    |     i     |    i+1    |
  // 	 |+         -|+         -|+         -|
  //  i-3/2       i-1/2       i+1/2       i+3/2
  //
  // ENO stencils:
  //
  //
  //                   |___________S2__________|
  //                   |                       |
  //           |___________S1__________|       |
  //           |                       |       |
  //   |___________S0__________|       |       |
  // ..|---o---|---o---|---o---|---o---|---o---|...
  //   | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
  //                                  -|
  //                                 i+1/2
  //
  //
  //                   |___________S0__________|
  //                   |                       |
  //                   |       |___________S1__________|
  //                   |       |                       |
  //                   |       |       |___________S2__________|
  //                 ..|---o---|---o---|---o---|---o---|---o---|...
  //                   | I{i-1}|  I{i} | I{i+1}| I{i+2}| I{i+3}|
  //                                   |+
  //                                 i+1/2
  //
  // WENO stencil: S{i} = [ I{i-2}, ..., I{i+3} ]
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  constexpr auto epsilon  = static_cast<sc_t>(1e-6);
  constexpr auto one	  = static_cast<sc_t>(1);
  constexpr auto two	  = static_cast<sc_t>(2);
  constexpr auto three	  = static_cast<sc_t>(3);
  constexpr auto four	  = two*two;
  constexpr auto five	  = three+two;
  constexpr auto six	  = three*two;
  constexpr auto seven	  = four+three;
  constexpr auto ten      = five*two;
  constexpr auto eleven   = five+six;
  constexpr auto twelve   = six*two;
  constexpr auto thirteen = six+seven;

  constexpr auto oneOvfour= one/four;
  constexpr auto oneOvsix = one/six;
  constexpr auto oneOvten = one/ten;
  constexpr auto threeOvten = three/ten;
  constexpr auto sixOvten = six/ten;
  constexpr auto thirteenOvtwelve = thirteen/twelve;

  {
    // u_{i+1/2}^{-}
    const auto p0 = (two*qim2 - seven*qim1 + eleven*qi)*oneOvsix;
    const auto p1 = ( -qim1   + five*qi    + two*qip1 )*oneOvsix;
    const auto p2 = (two*qi   + five*qip1  - qip2     )*oneOvsix;

    const auto B0 = thirteenOvtwelve* std::pow(qim2-two*qim1+qi  , two) +
		    oneOvfour*std::pow(qim2-four*qim1+three*qi, two);

    const auto B1 = thirteenOvtwelve* std::pow(qim1 -two*qi +qip1, two) +
                    oneOvfour*std::pow(qim1-qip1, two);

    const auto B2 = thirteenOvtwelve* std::pow(qi  -two*qip1+qip2, two) +
		    oneOvfour*std::pow(three*qi-four*qip1+qip2, two);

    const auto alpha0    = oneOvten/std::pow(epsilon   + B0, two);
    const auto alpha1    = sixOvten/std::pow(epsilon   + B1, two);
    const auto alpha2    = threeOvten/std::pow(epsilon + B2, two);
    const auto alphaSInv = one/(alpha0 + alpha1 + alpha2);

    const auto w0 = alpha0*alphaSInv;
    const auto w1 = alpha1*alphaSInv;
    const auto w2 = alpha2*alphaSInv;

    uNeg = w0*p0 + w1*p1 + w2*p2;
  }

  {
    // u_{i+1/2}^{+}
    const auto p0 = (-qim1        + five*qi    + two*qip1 )*oneOvsix;
    const auto p1 = (two*qi       + five*qip1  - qip2     )*oneOvsix;
    const auto p2 = (eleven*qip1  - seven*qip2 + two*qip3 )*oneOvsix;

    const auto B0 = thirteenOvtwelve* std::pow(qim1-two*qi+qip1   , two) +
		    oneOvfour*std::pow(qim1-four*qi+three*qip1, two);

    const auto B1 = thirteenOvtwelve* std::pow(qi  -two*qip1 +qip2, two) +
		    oneOvfour*std::pow(qi-qip2, two);

    const auto B2 = thirteenOvtwelve* std::pow(qip1-two*qip2+qip3 , two) +
		    oneOvfour*std::pow(three*qip1-four*qip2+qip3, two);

    const auto alpha0    = threeOvten/std::pow(epsilon + B0, two);
    const auto alpha1    = sixOvten/std::pow(epsilon   + B1, two);
    const auto alpha2    = oneOvten/std::pow(epsilon   + B2, two);
    const auto alphaSInv = one/(alpha0 + alpha1 + alpha2);

    const auto w0 = alpha0*alphaSInv;
    const auto w1 = alpha1*alphaSInv;
    const auto w2 = alpha2*alphaSInv;

    uPos = w0*p0 + w1*p1 + w2*p2;
  }
}

}}
#endif
