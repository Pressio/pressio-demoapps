
#ifndef PRESSIODEMOAPPS_ADV_DIFF_2D_RUSANOV_FLUX_VALUES_FUNC_HPP_
#define PRESSIODEMOAPPS_ADV_DIFF_2D_RUSANOV_FLUX_VALUES_FUNC_HPP_

namespace pressiodemoapps{ namespace impladvdiff2d{

template<class sc_t, typename T, class normal_t>
void burgers_rusanov_flux_2d(T & F,
			     const T & qL,
			     const T & qR,
			     const normal_t & n)
{
  constexpr sc_t one  = static_cast<sc_t>(1);
  constexpr sc_t two  = static_cast<sc_t>(2);
  constexpr sc_t half = one/two;

  const sc_t alpha_0 = std::max( std::abs(qL[0]), std::abs(qR[0]));
  const sc_t alpha_1 = std::max( std::abs(qL[1]), std::abs(qR[1]));

  F[0] = alpha_0*(qL[0] - qR[0]);
  F[0] += n[0]*(qL[0]*qL[0] + qR[0]*qR[0]);
  F[0] += n[1]*(qL[0]*qL[1] + qR[0]*qR[1]);
  F[0] *= half;

  F[1] = alpha_1*(qL[1] - qR[1]);
  F[1] += n[0]*(qL[0]*qL[1] + qR[0]*qR[1]);
  F[1] += n[1]*(qL[1]*qL[1] + qR[1]*qR[1]);
  F[1] *= half;
}

template<class sc_t, typename T1, class T2, class normal_t>
void burgers_rusanov_flux_jacobian_2d(T1 & JL, T1 & JR,
				      const T2 & qL,
				      const T2 & qR,
				      const normal_t & n)
{

  constexpr sc_t one  = static_cast<sc_t>(1);
  constexpr sc_t two  = static_cast<sc_t>(2);
  constexpr sc_t half = one/two;

  if (std::abs(qL[0]) > std::abs(qR[0])) {
    JL(0,0) = (2*qL[0]-qR[0])*std::copysign(1,qL[0]) + n[0]*2*qL[0] + n[1]*qL[1];
    JR(0,0) = n[0]*2*qR[0] + n[1]*qR[1] - std::abs(qL[0]);
  } else {
    JL(0,0) = n[0]*2*qL[0] + n[1]*qL[1] + std::abs(qR[0]);
    JR(0,0) = (qL[0]-2*qR[0])*std::copysign(1,qR[0]) + n[0]*2*qR[0] + n[1]*qR[1];
  }
  JL(0,0) *= 0.5;
  JR(0,0) *= 0.5;

  if (std::abs(qL[1]) > std::abs(qR[1])) {
    JL(1,1) = (2*qL[1]-qR[1])*std::copysign(1,qL[1]) + n[0]*qL[0] + n[1]*2*qL[1];
    JR(1,1) = n[0]*qR[0] + n[1]*2*qR[1] - std::abs(qL[1]);
  } else {
    JL(1,1) = n[0]*qL[0] + n[1]*2*qL[1] + std::abs(qR[1]);
    JR(1,1) = (qL[1]-2*qR[1])*std::copysign(1,qR[1]) + n[0]*qR[0] + n[1]*2*qR[1];
  }
  JL(1,1) *= 0.5;
  JR(1,1) *= 0.5;

  JL(0,1) = 0.5*n[1]*qL[0];
  JL(1,0) = 0.5*n[0]*qL[1];

  JR(0,1) = 0.5*n[1]*qR[0];
  JR(1,0) = 0.5*n[0]*qR[1];

}

}}
#endif
