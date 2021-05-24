
#ifndef PRESSIODEMOAPPS_WENO5IMPL_HPP_
#define PRESSIODEMOAPPS_WENO5IMPL_HPP_

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

  constexpr auto epsilon = 1e-6;

  auto vmm  = qim2;
  auto vm   = qim1;
  auto v    = qi;
  auto vp   = qip1;
  auto vpp  = qip2;

  // u_{i+1/2}^{-}
  const auto p0n = (2.*vmm - 7.*vm + 11.*v)/6.;
  const auto p1n = ( -vm   + 5.*v  + 2.*vp)/6.;
  const auto p2n = (2.*v   + 5.*vp - vpp )/6.;

  const auto B0n = (13./12.)* std::pow(vmm-2.*vm+v,  2.) + (1./4.)*std::pow(vmm-4.*vm+3.*v, 2.);
  const auto B1n = (13./12.)* std::pow(vm -2.*v +vp, 2.) + (1./4.)*std::pow(vm-vp, 2.);
  const auto B2n = (13./12.)* std::pow(v  -2.*vp+vpp,2.) + (1./4.)*std::pow(3.*v-4.*vp+vpp, 2.);

  constexpr auto d0n = 1./10.;
  constexpr auto d1n = 6./10.;
  constexpr auto d2n = 3./10.;
  const auto alpha0n = d0n/std::pow(epsilon + B0n, 2.);
  const auto alpha1n = d1n/std::pow(epsilon + B1n, 2.);
  const auto alpha2n = d2n/std::pow(epsilon + B2n, 2.);
  const auto alphasumn = alpha0n + alpha1n + alpha2n;

  const auto w0n = alpha0n/alphasumn;
  const auto w1n = alpha1n/alphasumn;
  const auto w2n = alpha2n/alphasumn;

  uNeg = w0n*p0n + w1n*p1n + w2n*p2n;

  // u_{i+1/2}^{+}
  vmm  = qim1;
  vm   = qi;
  v    = qip1;
  vp   = qip2;
  vpp  = qip3;

  const auto p0p = (-vmm   + 5.*vm + 2.*v   )/6.;
  const auto p1p = (2*vm   + 5.*v  - vp     )/6.;
  const auto p2p = (11.*v  - 7.*vp + 2.*vpp )/6.;

  const auto B0p = (13./12.)* std::pow(vmm-2.*vm+v,  2.) + (1./4.)*std::pow(vmm-4.*vm+3.*v, 2.);
  const auto B1p = (13./12.)* std::pow(vm -2.*v +vp, 2.) + (1./4.)*std::pow(vm-vp, 2.);
  const auto B2p = (13./12.)* std::pow(v  -2.*vp+vpp,2.) + (1./4.)*std::pow(3.*v-4.*vp+vpp, 2.);

  constexpr auto d0p = 3./10.;
  constexpr auto d1p = 6./10.;
  constexpr auto d2p = 1./10.;
  const auto alpha0p = d0p/std::pow(epsilon + B0p, 2.);
  const auto alpha1p = d1p/std::pow(epsilon + B1p, 2.);
  const auto alpha2p = d2p/std::pow(epsilon + B2p, 2.);
  const auto alphasump = alpha0p + alpha1p + alpha2p;

  const auto w0p = alpha0p/alphasump;
  const auto w1p = alpha1p/alphasump;
  const auto w2p = alpha2p/alphasump;

  uPos = w0p*p0p + w1p*p1p + w2p*p2p;
}

}}
#endif
