
#ifndef PRESSIODEMOAPPS_SWE_JACOBIANS_HPP_
#define PRESSIODEMOAPPS_SWE_JACOBIANS_HPP_

namespace pressiodemoapps{ namespace implswe{

template<class T, typename T2, typename sc_t, typename normal_t>
void swe_rusanov_flux_jacobian_three_dof(T & JL,
				    T & JR,
				    const T2 & qL,
				    const T2 & qR,
				    const normal_t & n,
				    const sc_t g)
{

// Computes the flux Jacobian for the shallow water equations
// INPUTS:
//    qL: conservative state vector in left cell
//    qR: conservative state vector in right cell
//    n: normal pointing from the left cell to the right cell

// OUTPUTS:
//  F   : the flux out of the left cell (into the right cell)

  const sc_t es = 1.e-30;

  const sc_t hL = qL(0);
  const sc_t uL = qL(1)/(hL + es);
  const sc_t vL = qL(2)/(hL + es);
  const sc_t unL = uL*n[0] + vL*n[1];

  const sc_t hR = qR(0);
  const sc_t uR = qR(1)/(hR + es);
  const sc_t vR = qR(2)/(hR + es);
  const sc_t unR = uR*n[0] + vR*n[1];

  // rho average
  const sc_t hm = 0.5*(hL + hR);
  const sc_t um = (unL*std::pow(hL,0.5) + unR*std::pow(hR,0.5) )/( std::pow(hL,0.5) + std::pow(hR,0.5)  + es);
  // eigenvalues
  const sc_t smax = abs(um) + abs(std::pow(g*hm, 0.5));

  const sc_t termL = (n[0]*qL(1) + n[1]*qL(2)) / std::pow(qL(0),2.);
  const sc_t termR = (n[0]*qR(1) + n[1]*qR(2)) / std::pow(qR(0),2.);
  const sc_t hL_sqrt = std::pow(hL,0.5);
  const sc_t hR_sqrt = std::pow(hR,0.5);

  const sc_t hsqrt_un = hL_sqrt*unL + hR_sqrt*unR + es;
  sc_t dsmaxL[3];
  sc_t dsmaxR[3];
  dsmaxL[0] = - std::abs( hsqrt_un ) / (2.*hL_sqrt*std::pow(hL_sqrt + hR_sqrt,2.) ) +
              (0.5*unL / hL_sqrt - hL_sqrt*termL )*hsqrt_un  / ( (hL_sqrt + hR_sqrt)*std::abs( hsqrt_un  ) ) +
              g/(std::pow(2.,3./2.)*std::pow(g*(hL + hR),0.5) );
  dsmaxL[1] = n[0]*hsqrt_un / ( hL_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );
  dsmaxL[2] = n[1]*hsqrt_un / ( hL_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );


  dsmaxR[0] = - std::abs( hsqrt_un ) / (2.*hR_sqrt*std::pow(hL_sqrt + hR_sqrt,2.) ) +
              (0.5*unR / hR_sqrt - hR_sqrt*termR )*hsqrt_un  / ( (hL_sqrt + hR_sqrt)*std::abs( hsqrt_un ) ) +
              g/(std::pow(2.,3./2.)*std::pow(g*(hL + hR),0.5) );
  dsmaxR[1] = n[0]*hsqrt_un / ( hR_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );
  dsmaxR[2] = n[1]*hsqrt_un / ( hR_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );

  // jacobian w.r.p to the left state
  JL(0,0) = -0.5*dsmaxL[0]*(qR(0) - qL(0)) + 0.5*(n[0]*uL + n[1]*vL - qL(0)*termL)  + 0.5*smax;
  JL(0,1) = 0.5*n[0] - 0.5*dsmaxL[1]*(qR(0) - qL(0));
  JL(0,2) = 0.5*n[1] - 0.5*dsmaxL[2]*(qR(0) - qL(0));

  JL(1,0) = 0.5*(g*n[0]*qL(0) - qL(1)*termL) - 0.5*dsmaxL[0]*(qR(1) - qL(1));
  JL(1,1) = n[0]*uL  + 0.5*n[1]*vL + 0.5*smax -0.5*dsmaxL[1]*(qR(1) - qL(1));
  JL(1,2) = 0.5*n[1]*uL  -0.5*dsmaxL[2]*(qR(1) - qL(1));

  JL(2,0) = 0.5*(g*n[1]*qL(0) - qL(2)*termL) - 0.5*dsmaxL[0]*(qR(2) - qL(2));
  JL(2,1) = 0.5*n[0]*vL  -0.5*dsmaxL[1]*(qR(2) - qL(2));
  JL(2,2) = n[1]*vL + 0.5*n[0]*uL + 0.5*smax -0.5*dsmaxL[2]*(qR(2) - qL(2));

   // jacobian w.r.p to the right state
  JR(0,0) = -0.5*dsmaxR[0]*(qR(0) - qL(0)) + 0.5*(n[0]*uR + n[1]*vR - qR(0)*termR)  - 0.5*smax;
  JR(0,1) = 0.5*n[0] - 0.5*dsmaxR[1]*(qR(0) - qL(0));
  JR(0,2) = 0.5*n[1] - 0.5*dsmaxR[2]*(qR(0) - qL(0));

  JR(1,0) = 0.5*(g*n[0]*qR(0) - qR(1)*termR) - 0.5*dsmaxR[0]*(qR(1) - qL(1));
  JR(1,1) = n[0]*uR + 0.5*n[1]*vR - 0.5*smax - 0.5*dsmaxR[1]*(qR(1) - qL(1));
  JR(1,2) = 0.5*n[1]*uR -0.5*dsmaxR[2]*(qR(1) - qL(1));

  JR(2,0) = 0.5*(g*n[1]*qR(0)  - qR(2)*termR) - 0.5*dsmaxR[0]*(qR(2) - qL(2));
  JR(2,1) = 0.5*n[0]*vR - 0.5*dsmaxR[1]*(qR(2) - qL(2));
  JR(2,2) = n[1]*vR + 0.5*n[0]*uR - 0.5*smax - 0.5*dsmaxR[2]*(qR(2) - qL(2));

}
}}
#endif
