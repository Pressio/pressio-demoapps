
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

  const sc_t epsilon  = static_cast<sc_t>(1.e-6);
  const sc_t one	  = static_cast<sc_t>(1.);
  const sc_t two	  = static_cast<sc_t>(2.);
  const sc_t three	  = static_cast<sc_t>(3.);
  const sc_t oneOvtwo   = one/two;
  const sc_t oneOvthree = one/three;
  const sc_t twoOvthree = two/three;

  {
    // u_{i+1/2}^{-}
    const auto p0 = (-qim1 + three*qi)*oneOvtwo;
    const auto p1 = (qi   + qip1)*oneOvtwo;

    const auto B0 = (qim1-qi)*(qim1 - qi);
    const auto B1 = (qi-qip1)*(qi-qip1);

    const auto alpha0    = oneOvthree/(epsilon*epsilon + two*epsilon*B0 + B0*B0);
    const auto alpha1    = twoOvthree/(epsilon*epsilon + two*epsilon*B1 + B1*B1);
    const auto alphaSInv = one/(alpha0 + alpha1);

    const auto w0 = alpha0*alphaSInv;
    const auto w1 = alpha1*alphaSInv;
    uNeg = w0*p0 + w1*p1;
  }

  {
    // u_{i+1/2}^{+}
    const auto p0 = ( qi + qip1)*oneOvtwo;
    const auto p1 = ( three*qip1 - qip2)*oneOvtwo;

    const auto B0 = (qi-qip1)*(qi-qip1);
    const auto B1 = (qip1-qip2)*(qip1-qip2);

    const auto alpha0    = twoOvthree/(epsilon*epsilon + two*epsilon*B0 + B0*B0);
    const auto alpha1    = oneOvthree/(epsilon*epsilon + two*epsilon*B1 + B1*B1);
    const auto alphaSInv = one/(alpha0 + alpha1);

    const auto w0 = alpha0*alphaSInv;
    const auto w1 = alpha1*alphaSInv;
    uPos = w0*p0 + w1*p1;
  }

}

template<typename sc_t, typename grad_t>
void weno3WithGrad(sc_t & uNeg,
		   sc_t & uPos,
		   grad_t & duNeg_dq,
		   grad_t & duPos_dq,
		   const sc_t & qim1,
		   const sc_t & qi,
		   const sc_t & qip1,
		   const sc_t & qip2)
{

  const sc_t epsilon  = static_cast<sc_t>(1e-6);
  const sc_t one	  = static_cast<sc_t>(1);
  const sc_t two	  = static_cast<sc_t>(2);
  const sc_t three	  = static_cast<sc_t>(3);
  const sc_t oneOvtwo   = one/two;
  const sc_t oneOvthree = one/three;
  const sc_t twoOvthree = two/three;

  {
    // u_{i+1/2}^{-}
    sc_t dp0_dq[3];
    sc_t dp1_dq[3];
    sc_t dalpha0_dq[3];
    sc_t dalpha1_dq[3];
    sc_t dalphaSInv_dq[3];
    const auto p0 = (-qim1 + three*qi)*oneOvtwo;
    const auto p1 = (qi   + qip1)*oneOvtwo;

    dp0_dq[0] = -1./2.;
    dp0_dq[1] = 3./2.;
    dp0_dq[2] = 0.;

    dp1_dq[0] = 0.;
    dp1_dq[1] = 1./2.;
    dp1_dq[2] = 1./2.;

    const auto B0 = (qim1-qi)*(qim1-qi);
    const auto B1 = (qi-qip1)*(qi-qip1);
    const auto alpha0    = oneOvthree/(epsilon*epsilon + two*epsilon*B0 + B0*B0);
    const auto alpha1    = twoOvthree/(epsilon*epsilon + two*epsilon*B1 + B1*B1);
    const auto alphaSInv = one/(alpha0 + alpha1);

    dalpha0_dq[0] = -4.*(qim1 - qi)/(3. * std::pow( std::pow(qim1 - qi,2.) + epsilon,3));
    dalpha0_dq[1] = 4.*(qim1 - qi)/(3.*std::pow(std::pow(qim1-qi,2.) + epsilon,3.));
    dalpha0_dq[2] = 0.;

    dalpha1_dq[0] = 0.;
    dalpha1_dq[1] = -8.*(qi - qip1)/(3.*std::pow(std::pow(qi-qip1,2.) + epsilon , 3.));
    dalpha1_dq[2] = 8.*(qi - qip1)/(3.*std::pow(std::pow(qi-qip1,2.) + epsilon , 3.));

    dalphaSInv_dq[0] = (4.*(qim1-qi))/(3.*std::pow(std::pow(qim1-qi,2.)+epsilon,3.)*std::pow(2./(3.*std::pow(std::pow(qi-qip1,2.)+epsilon,2.))+1./(3.*std::pow(std::pow(qim1-qi,2.)+epsilon,2.)),2.));
    dalphaSInv_dq[1] = -((4.*(qim1-qi))/(3.*std::pow(std::pow(qim1-qi,2.)+epsilon,3.))-(8.*(qi-qip1))/(3.*std::pow(std::pow(qi-qip1,2.)+epsilon,3.)))/std::pow(2./(3.*std::pow(std::pow(qi-qip1,2.)+epsilon,2.))+1./(3.*std::pow(std::pow(qim1-qi,2.)+epsilon,2.)),2.);
    dalphaSInv_dq[2] = -(8.*(qi-qip1))/(3.*std::pow(2./(3.*std::pow(std::pow(qi-qip1,2.)+epsilon,2.))+1./(3.*std::pow(std::pow(qim1-qi,2.)+epsilon,2.)),2.)*std::pow(std::pow(qi-qip1,2.)+epsilon,3.));
    const auto w0 = alpha0*alphaSInv;
    const auto w1 = alpha1*alphaSInv;

    for (int i = 0; i < 3; i++){
      duNeg_dq[i] = (dalpha0_dq[i]*alphaSInv + dalphaSInv_dq[i]*alpha0)*p0 + dp0_dq[i]*w0;
      duNeg_dq[i] += (dalpha1_dq[i]*alphaSInv + dalphaSInv_dq[i]*alpha1)*p1 + dp1_dq[i]*w1;
    }
    duNeg_dq[3] = 0.;
    uNeg = w0*p0 + w1*p1;

  }
  {
    // u_{i+1/2}^{+}
    sc_t dp0_dq[3];
    sc_t dp1_dq[3];
    sc_t dalpha0_dq[3];
    sc_t dalpha1_dq[3];
    sc_t dalphaSInv_dq[3];

    const auto p0 = ( qi + qip1)*oneOvtwo;
    const auto p1 = ( three*qip1 - qip2)*oneOvtwo;
    const auto B0 = (qi-qip1)*(qi-qip1);
    const auto B1 = (qip1-qip2)*(qip1-qip2);
    const auto alpha0    = twoOvthree/(epsilon*epsilon + two*epsilon*B0 + B0*B0);
    const auto alpha1    = oneOvthree/(epsilon*epsilon + two*epsilon*B1 + B1*B1);
    const auto alphaSInv = one/(alpha0 + alpha1);
    const auto w0 = alpha0*alphaSInv;
    const auto w1 = alpha1*alphaSInv;

    dp0_dq[0] = 1./2.;
    dp0_dq[1] = 1/2.;
    dp0_dq[2] = 0.;

    dp1_dq[0] = 0.;
    dp1_dq[1] = 3./2.;
    dp1_dq[2] = -1./2.;

    dalpha0_dq[0] = -8.*(qi - qip1)/(3.*std::pow(std::pow(qi-qip1,2.)+epsilon,3.));
    dalpha0_dq[1] = 8.*(qi - qip1)/(3.*std::pow(std::pow(qi-qip1,2.)+epsilon,3.));
    dalpha0_dq[2] = 0.;

    dalpha1_dq[0] = 0.;
    dalpha1_dq[1] = -4.*(qip1 - qip2)/(3.*std::pow(std::pow(qip1-qip2,2.) + epsilon,3.));
    dalpha1_dq[2] = 4.*(qip1 - qip2)/(3. * std::pow(std::pow(qip1 - qip2,2)+epsilon,3.));

    dalphaSInv_dq[0] = 8.*(qi - qip1)/(3.*std::pow(std::pow(qi-qip1,2.) + epsilon,3.)*std::pow(1./(3.*std::pow(std::pow(qip1-qip2,2.) + epsilon,2.))+2./(3.*std::pow(std::pow(qi-qip1,2.) + epsilon,2.)),2.));
    dalphaSInv_dq[1] = -(8.*(qi - qip1)/(3*std::pow(std::pow(qi-qip1,2.) + epsilon,3.))-4.*(qip1 - qip2)/(3.*std::pow(std::pow(qip1-qip2,2.) + epsilon,3.)))/std::pow(1./(3.*std::pow(std::pow(qip1-qip2,2.) + epsilon,2.))+2./(3.*std::pow(std::pow(qi-qip1,2.) + epsilon,2.)),2.);
    dalphaSInv_dq[2] = -4.*(qip1 - qip2)/(3.*std::pow(1./(3.*std::pow(std::pow(qip1-qip2,2.) + epsilon,2.))+2./(3.*std::pow(std::pow(qi-qip1,2.) + epsilon,2.)),2.)*std::pow(std::pow(qip1-qip2,2.) + epsilon,3.));

    for (int i = 0; i < 3; i++){
      duPos_dq[i+1] = (dalpha0_dq[i]*alphaSInv + dalphaSInv_dq[i]*alpha0)*p0 + dp0_dq[i]*w0;
      duPos_dq[i+1] += (dalpha1_dq[i]*alphaSInv + dalphaSInv_dq[i]*alpha1)*p1 + dp1_dq[i]*w1;
    }
    duPos_dq[0] = 0.;
    uPos = w0*p0 + w1*p1;

  }

}



}}
#endif
