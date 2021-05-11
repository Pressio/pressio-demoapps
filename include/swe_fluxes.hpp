
#ifndef ODSWE_RUSANOV_FLUX_HPP_
#define ODSWE_RUSANOV_FLUX_HPP_

#include <math.h>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <iostream>
#include <fstream>

/*
F \in R^3, flux
q \in R^N, state
isL \in NaturalNumbers, left index
isR \in NaturalNumbers, right index
n \in R^2, normalv ector
g \in R, gravity
*/
template<typename flux_t,typename state_t, typename normal_t, typename scalar_t>
void sweRusanovFluxFullStateIn(flux_t & F,
			       const state_t & q,
			       const int isL,
			       const int isR,
			       const normal_t & n,
			       const scalar_t g)
{
  // Computes the flux for the shallow water equations
  const scalar_t es = 1.e-30;

  // left
  const auto hL = q(isL + 0);
  const auto uL = q(isL + 1)/(hL + es);
  const auto vL = q(isL + 2)/(hL + es);
  const auto unL = uL*n[0] + vL*n[1];
  const auto pL = 0.5*g*std::pow(hL,2.);
  scalar_t FL[ 3 ];
  FL[0] = hL*unL;
  FL[1] = q(isL + 1)*unL + pL*n[0];
  FL[2] = q(isL + 2)*unL + pL*n[1];

  // right
  const auto hR = q(isR + 0);
  const auto uR = q(isR + 1)/(hR + es);
  const auto vR = q(isR + 2)/(hR + es);
  const auto unR = uR*n[0] + vR*n[1];
  const auto pR = 0.5*g*std::pow(hR,2.);
  scalar_t FR[3];
  FR[0] = hR*unR;
  FR[1] = q(isR + 1)*unR + pR*n[0];
  FR[2] = q(isR + 2)*unR + pR*n[1];

  const auto hm = 0.5*(hL + hR);
  const auto umNum = unL*std::pow(hL,0.5) + unR*std::pow(hR,0.5);
  const auto umDen = ( std::pow(hL,0.5) + std::pow(hR,0.5)  + es );
  const auto um = umNum/umDen;
  const auto smax = abs(um) + abs(std::pow(g*hm, 0.5));

  // flux assembly
  F[0] = 0.5*(FL[0]+FR[0])-0.5*smax*(q(isR + 0) - q(isL + 0));
  F[1] = 0.5*(FL[1]+FR[1])-0.5*smax*(q(isR + 1) - q(isL + 1));
  F[2] = 0.5*(FL[2]+FR[2])-0.5*smax*(q(isR + 2) - q(isL + 2));
}


/* Rusanov flux Jacobian that takes the entire state vector as an input
   JL \in R^{3x3}, left jacobian
   JR \in R^{3x3}, right jacobian
   q \in R^N, state
   isL \in NaturalNumbers, left index
   isR \in NaturalNumbers, right index
   n \in R^2, normalv ector
   g \in R, gravity
*/
template<typename jac_t,typename state_t, typename normal_t, typename scalar_t>
void rusanovFluxJacobianFullStateIn(jac_t & JL,
				    jac_t & JR,
				    const state_t & q,
				    const int gidL,
				    const int gidR,
				    const normal_t & n,
				    const scalar_t g)
{
  const scalar_t es = 1.e-30;
  const scalar_t hL = q(gidL + 0);
  const scalar_t uL = q(gidL + 1)/(hL + es);
  const scalar_t vL = q(gidL + 2)/(hL + es);
  const scalar_t unL = uL*n[0] + vL*n[1];

  const scalar_t hR = q(gidR + 0);
  const scalar_t uR = q(gidR + 1)/(hR + es);
  const scalar_t vR = q(gidR + 2)/(hR + es);
  const scalar_t unR = uR*n[0] + vR*n[1];

  // rho average
  const scalar_t hm = 0.5*(hL + hR);
  const scalar_t um = (unL*std::pow(hL,0.5) + unR*std::pow(hR,0.5) )/( std::pow(hL,0.5) + std::pow(hR,0.5)  + es);
  // eigenvalues
  const scalar_t smax = abs(um) + abs(std::pow(g*hm, 0.5));

  const scalar_t termL = (n[0]*q(gidL + 1) + n[1]*q(gidL + 2)) / std::pow(q(gidL + 0),2.);
  const scalar_t termR = (n[0]*q(gidR + 1) + n[1]*q(gidR + 2)) / std::pow(q(gidR + 0),2.);
  const scalar_t hL_sqrt = std::pow(hL,0.5);
  const scalar_t hR_sqrt = std::pow(hR,0.5);

  const scalar_t hsqrt_un = hL_sqrt*unL + hR_sqrt*unR + es;
  scalar_t dsmaxL[3];
  scalar_t dsmaxR[3];
  dsmaxL[0] = - std::abs(hsqrt_un) / (2.*hL_sqrt*std::pow(hL_sqrt + hR_sqrt,2.) ) +
    (0.5*unL / hL_sqrt - hL_sqrt*termL )*hsqrt_un  / ( (hL_sqrt + hR_sqrt)*std::abs( hsqrt_un  ) ) +
    g/(std::pow(2.,3./2.)*std::pow(g*(hL + hR),0.5) );
  dsmaxL[1] = n[0]*hsqrt_un / ( hL_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );
  dsmaxL[2] = n[1]*hsqrt_un / ( hL_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );

  dsmaxR[0] = - std::abs(hsqrt_un) / (2.*hR_sqrt*std::pow(hL_sqrt + hR_sqrt,2.) ) +
    (0.5*unR / hR_sqrt - hR_sqrt*termR )*hsqrt_un  / ( (hL_sqrt + hR_sqrt)*std::abs( hsqrt_un ) ) +
    g/(std::pow(2.,3./2.)*std::pow(g*(hL + hR),0.5) );
  dsmaxR[1] = n[0]*hsqrt_un / ( hR_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );
  dsmaxR[2] = n[1]*hsqrt_un / ( hR_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );

  // jacobian w.r.p to the left state
  JL[0][0] = -0.5*dsmaxL[0]*(q(gidR + 0) - q(gidL + 0)) + 0.5*(n[0]*uL + n[1]*vL - q(gidL + 0)*termL)  + 0.5*smax;
  JL[0][1] = 0.5*n[0] - 0.5*dsmaxL[1]*(q(gidR + 0) - q(gidL + 0));
  JL[0][2] = 0.5*n[1] - 0.5*dsmaxL[2]*(q(gidR + 0) - q(gidL + 0));

  JL[1][0] = 0.5*(g*n[0]*q(gidL + 0) - q(gidL + 1)*termL) - 0.5*dsmaxL[0]*(q(gidR + 1) - q(gidL + 1));
  JL[1][1] = n[0]*uL  + 0.5*n[1]*vL + 0.5*smax -0.5*dsmaxL[1]*(q(gidR + 1) - q(gidL + 1));
  JL[1][2] = 0.5*n[1]*uL  -0.5*dsmaxL[2]*(q(gidR + 1) - q(gidL + 1));

  JL[2][0] = 0.5*(g*n[1]*q(gidL + 0) - q(gidL + 2)*termL) - 0.5*dsmaxL[0]*(q(gidR + 2) - q(gidL + 2));
  JL[2][1] = 0.5*n[0]*vL  -0.5*dsmaxL[1]*(q(gidR + 2) - q(gidL + 2));
  JL[2][2] = n[1]*vL + 0.5*n[0]*uL + 0.5*smax -0.5*dsmaxL[2]*(q(gidR + 2) - q(gidL + 2));

  // jacobian w.r.p to the right state
  JR[0][0] = -0.5*dsmaxR[0]*(q(gidR + 0) - q(gidL + 0)) + 0.5*(n[0]*uR + n[1]*vR - q(gidR + 0)*termR)  - 0.5*smax;
  JR[0][1] = 0.5*n[0] - 0.5*dsmaxR[1]*(q(gidR + 0) - q(gidL + 0));
  JR[0][2] = 0.5*n[1] - 0.5*dsmaxR[2]*(q(gidR + 0) - q(gidL + 0));

  JR[1][0] = 0.5*(g*n[0]*q(gidR + 0) - q(gidR + 1)*termR) - 0.5*dsmaxR[0]*(q(gidR + 1) - q(gidL + 1));
  JR[1][1] = n[0]*uR + 0.5*n[1]*vR - 0.5*smax - 0.5*dsmaxR[1]*(q(gidR + 1) - q(gidL + 1));
  JR[1][2] = 0.5*n[1]*uR -0.5*dsmaxR[2]*(q(gidR + 1) - q(gidL + 1));

  JR[2][0] = 0.5*(g*n[1]*q(gidR + 0)  - q(gidR + 2)*termR) - 0.5*dsmaxR[0]*(q(gidR + 2) - q(gidL + 2));
  JR[2][1] = 0.5*n[0]*vR - 0.5*dsmaxR[1]*(q(gidR + 2) - q(gidL + 2));
  JR[2][2] = n[1]*vR + 0.5*n[0]*uR - 0.5*smax - 0.5*dsmaxR[2]*(q(gidR + 2) - q(gidL + 2));

}


template<typename sc_t, class normal_t>
void sweLF(std::array<sc_t,3> & F,
	const std::array<sc_t,3> & qL,
	const std::array<sc_t,3> & qR,
	const normal_t & n,
	const sc_t g)
{
  const sc_t es = 1.e-30;

  // left
  const auto hL = qL[0];
  const auto uL = qL[1]/(hL + es);
  const auto vL = qL[2]/(hL + es);
  const auto unL = uL*n[0] + vL*n[1];
  const auto pL = 0.5*g*std::pow(hL,2.);
  sc_t FL[ 3 ];
  FL[0] = hL*unL;
  FL[1] = qL[1]*unL + pL*n[0];
  FL[2] = qL[2]*unL + pL*n[1];

  // right
  const auto hR = qR[0];
  const auto uR = qR[1]/(hR + es);
  const auto vR = qR[2]/(hR + es);
  const auto unR = uR*n[0] + vR*n[1];
  const auto pR = 0.5*g*std::pow(hR,2.);
  sc_t FR[3];
  FR[0] = hR*unR;
  FR[1] = qR[1]*unR + pR*n[0];
  FR[2] = qR[2]*unR + pR*n[1];

  // rho average
  const auto hm = 0.5*(hL + hR);
  const auto umNum = unL*std::pow(hL,0.5) + unR*std::pow(hR,0.5);
  const auto umDen =    (std::pow(hL,0.5) + std::pow(hR,0.5) + es);
  const auto um = umNum/umDen;
  // eigenvalues
  const auto smax = abs(um) + abs(std::pow(g*hm, 0.5));
  // flux assembly
  F[0] = 0.5*(FL[0]+FR[0])-0.5*smax*(qR[0] - qL[0]);
  F[1] = 0.5*(FL[1]+FR[1])-0.5*smax*(qR[1] - qL[1]);
  F[2] = 0.5*(FL[2]+FR[2])-0.5*smax*(qR[2] - qL[2]);
}


template<typename sc_t, typename state_t>
void weno5(std::array<sc_t,3> & uNeg,
	   std::array<sc_t,3> & uPos,
	   const state_t & q,
	   const int im2,
	   const int im1,
	   const int i,
	   const int ip1,
	   const int ip2,
	   const int ip3)
{
  // Domain cells (I{i}) reference:
  //
  //                |           |   u(i)    |           |
  //                |  u(i-1)   |___________|           |
  //                |___________|           |   u(i+1)  |
  //                |           |           |___________|
  //             ...|-----0-----|-----0-----|-----0-----|...
  //                |    i-1    |     i     |    i+1    |
  //                |+         -|+         -|+         -|
  //              i-3/2       i-1/2       i+1/2       i+3/2
  //
  // ENO stencils (S{r}) reference:
  //
  //
  //                           |___________S2__________|
  //                           |                       |
  //                   |___________S1__________|       |
  //                   |                       |       |
  //           |___________S0__________|       |       |
  //         ..|---o---|---o---|---o---|---o---|---o---|...
  //           | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
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
  // WENO stencil: S{i} = [ I{i-2},...,I{i+3} ]
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  const auto epsilon = 1e-6;

  for (int dof=0; dof<3; ++dof)
  {
    auto vmm  = q[im2+dof];
    auto vm   = q[im1+dof];
    auto v    = q[i+dof];
    auto vp   = q[ip1+dof];
    auto vpp  = q[ip2+dof];

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

    uNeg[dof] = w0n*p0n + w1n*p1n + w2n*p2n;

    // u_{i+1/2}^{+}
    vmm  = q[im1+dof];
    vm   = q[i+dof];
    v    = q[ip1+dof];
    vp   = q[ip2+dof];
    vpp  = q[ip3+dof];

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

    uPos[dof] = w0p*p0p + w1p*p1p + w2p*p2p;
  }
}


// template<typename flux_t,typename state_t, typename normal_t, typename scalar_t>
// void rusanovFluxFullStateInSplitDofs( flux_t  & F,
// 				      const state_t & q,
// 				      const int isL,
// 				      const int isR,
// 				      const normal_t n,
// 				      const scalar_t g)
// {
//   const scalar_t es = 1.e-30;

//   // left
//   const auto hL = q(isL, 0);
//   const auto uL = q(isL, 1)/(hL + es);
//   const auto vL = q(isL, 2)/(hL + es);
//   const auto unL = uL*n[0] + vL*n[1];
//   const auto pL = 0.5*g*std::pow(hL,2.);
//   scalar_t FL[ 3 ];
//   FL[0] = hL*unL;
//   FL[1] = q(isL,1)*unL + pL*n[0];
//   FL[2] = q(isL,2)*unL + pL*n[1];

//   // right
//   const auto hR = q(isR,0);
//   const auto uR = q(isR,1)/(hR + es);
//   const auto vR = q(isR,2)/(hR + es);
//   const auto unR = uR*n[0] + vR*n[1];
//   const auto pR = 0.5*g*std::pow(hR,2.);
//   scalar_t FR[3];
//   FR[0] = hR*unR;
//   FR[1] = q(isR,1)*unR + pR*n[0];
//   FR[2] = q(isR,2)*unR + pR*n[1];

//   // rho average
//   const auto hm = 0.5*(hL + hR);
//   const auto umNum = unL*std::pow(hL,0.5) + unR*std::pow(hR,0.5);
//   const auto umDen = ( std::pow(hL,0.5) + std::pow(hR,0.5)  + es );
//   const auto um = umNum/umDen;
//   // eigenvalues
//   const auto smax = abs(um) + abs(std::pow(g*hm, 0.5));
//   // flux assembly
//   F[0] = 0.5*(FL[0]+FR[0])-0.5*smax*(q(isR,0) - q(isL,0));
//   F[1] = 0.5*(FL[1]+FR[1])-0.5*smax*(q(isR,1) - q(isL,1));
//   F[2] = 0.5*(FL[2]+FR[2])-0.5*smax*(q(isR,2) - q(isL,2));
// }

// /* Same as above, but takes in the left state and right state as inputs
// F \in R^3, flux
// qL \in R^3, left state
// qR \in R^3, right state
// n \in R^2, normalv ector
// g \in R, gravity
// */
// template<typename flux_t,typename state_t, typename normal_t, typename scalar_t>
// void rusanovFluxLeftRightStateIn( flux_t  & F,
// 				  const state_t & qL,
// 				  const state_t & qR,
// 				  const normal_t n,
// 				  const scalar_t g)
// {
// // Computes the flux for the shallow water equations
//   const scalar_t es = 1.e-30;
//   const scalar_t hL = qL[0];
//   const scalar_t uL = qL[1]/(hL + es);
//   const scalar_t vL = qL[2]/(hL + es);
//   const scalar_t unL = uL*n[0] + vL*n[1];

//   const scalar_t pL = 0.5*g*std::pow(hL,2.);
//   scalar_t FL[ 3 ];
//   FL[0] = hL*unL;
//   FL[1] = qL[1]*unL + pL*n[0];
//   FL[2] = qL[2]*unL + pL*n[1];

//   const scalar_t hR = qR[0];
//   const scalar_t uR = qR[1]/(hR + es);
//   const scalar_t vR = qR[2]/(hR + es);
//   const scalar_t unR = uR*n[0] + vR*n[1];
//   const scalar_t pR = 0.5*g*std::pow(hR,2.);
//   scalar_t FR[3];
//   FR[0] = hR*unR;
//   FR[1] = qR[1]*unR + pR*n[0];
//   FR[2] = qR[2]*unR + pR*n[1];

//   // rho average
//   const scalar_t hm = 0.5*(hL + hR);
//   const scalar_t um = (unL*std::pow(hL,0.5) + unR*std::pow(hR,0.5) )/( std::pow(hL,0.5) + std::pow(hR,0.5)  + es);
//   // eigenvalues
//   const scalar_t smax = abs(um) + abs(std::pow(g*hm, 0.5));
//   // flux assembly
//   F[0]    = 0.5*(FL[0]+FR[0])-0.5*smax*(qR[0] - qL[0]);
//   F[1]    = 0.5*(FL[1]+FR[1])-0.5*smax*(qR[1] - qL[1]);
//   F[2]    = 0.5*(FL[2]+FR[2])-0.5*smax*(qR[2] - qL[2]);
// }

#endif
