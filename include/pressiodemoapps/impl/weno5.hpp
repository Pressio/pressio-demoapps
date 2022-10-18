/*
//@HEADER
// ************************************************************************
//
// weno5.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
  // 	  |           |   u(i)    |           |
  // 	  |  u(i-1)   |___________|           |
  // 	  |___________|           |   u(i+1)  |
  // 	  |           |           |___________|
  // ...|-----0-----|-----0-----|-----0-----|...
  // 	  |    i-1    |     i     |    i+1    |
  // 	  |+         -|+         -|+         -|
  //   i-3/2       i-1/2       i+1/2       i+3/2
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
  //                          -|
  //                          i+1/2
  //
  //
  //           |___________S0__________|
  //           |                       |
  //           |       |___________S1__________|
  //           |       |                       |
  //           |       |       |___________S2__________|
  //         ..|---o---|---o---|---o---|---o---|---o---|...
  //           | I{i-1}|  I{i} | I{i+1}| I{i+2}| I{i+3}|
  //                           |+
  //                          i+1/2
  //
  // WENO stencil: S{i} = [ I{i-2}, ..., I{i+3} ]
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  const auto epsilon  = static_cast<sc_t>(1e-6);
  const auto one	  = static_cast<sc_t>(1);
  const auto two	  = static_cast<sc_t>(2);
  const auto three	  = static_cast<sc_t>(3);
  const auto four	  = two*two;
  const auto five	  = three+two;
  const auto six	  = three*two;
  const auto seven	  = four+three;
  const auto ten      = five*two;
  const auto eleven   = five+six;
  const auto twelve   = six*two;
  const auto thirteen = six+seven;

  const auto oneOvfour= one/four;
  const auto oneOvsix = one/six;
  const auto oneOvten = one/ten;
  const auto threeOvten = three/ten;
  const auto sixOvten = six/ten;
  const auto thirteenOvtwelve = thirteen/twelve;

  {
    // u_{i+1/2}^{-}
    const auto p0 = (two*qim2 - seven*qim1 + eleven*qi)*oneOvsix;
    const auto p1 = ( -qim1   + five*qi    + two*qip1 )*oneOvsix;
    const auto p2 = (two*qi   + five*qip1  - qip2     )*oneOvsix;

    const auto B0 = thirteenOvtwelve* (qim2 - two*qim1 + qi)*(qim2 - two*qim1 + qi) +
		    oneOvfour*(qim2-four*qim1+three*qi)*(qim2-four*qim1+three*qi);

    const auto B1 = thirteenOvtwelve* (qim1 -two*qi +qip1)*(qim1 -two*qi +qip1) +
                    oneOvfour*(qim1-qip1)*(qim1-qip1);

    const auto B2 = thirteenOvtwelve* (qi  -two*qip1+qip2)*(qi  -two*qip1+qip2) +
		    oneOvfour*(three*qi-four*qip1+qip2)*(three*qi-four*qip1+qip2);

    const auto alpha0    = oneOvten/(epsilon*epsilon + 2.*epsilon*B0 + B0*B0);
    const auto alpha1    = sixOvten/(epsilon*epsilon + 2.*epsilon*B1 + B1*B1);
    const auto alpha2    = threeOvten/(epsilon*epsilon + 2.*epsilon*B2 + B2*B2);
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

    const auto B0 = thirteenOvtwelve*(qim1-two*qi+qip1)*(qim1-two*qi+qip1) +
		    oneOvfour*(qim1-four*qi+three*qip1)*(qim1-four*qi+three*qip1);

    const auto B1 = thirteenOvtwelve*(qi  -two*qip1 +qip2)*(qi  -two*qip1 +qip2)+
		    oneOvfour*(qi-qip2)*(qi-qip2);

    const auto B2 = thirteenOvtwelve*(qip1-two*qip2+qip3)*(qip1-two*qip2+qip3) +
		    oneOvfour*(three*qip1-four*qip2+qip3)*(three*qip1-four*qip2+qip3);

    const auto alpha0    = threeOvten/(epsilon*epsilon + 2.*epsilon*B0 + B0*B0);
    const auto alpha1    = sixOvten/(epsilon*epsilon + 2.*epsilon*B1 + B1*B1);
    const auto alpha2    = oneOvten/(epsilon*epsilon + 2.*epsilon*B2 + B2*B2);
    const auto alphaSInv = one/(alpha0 + alpha1 + alpha2);

    const auto w0 = alpha0*alphaSInv;
    const auto w1 = alpha1*alphaSInv;
    const auto w2 = alpha2*alphaSInv;

    uPos = w0*p0 + w1*p1 + w2*p2;
  }
}

template<typename sc_t, typename jac_t>
void weno5WithGrad(sc_t & uNeg,
     sc_t & uPos,
     jac_t & duNeg_dq,
	   jac_t & duPos_dq,
	   const sc_t & qim2,
	   const sc_t & qim1,
	   const sc_t & qi,
	   const sc_t & qip1,
	   const sc_t & qip2,
	   const sc_t & qip3)
{

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  // 	  |           |   u(i)    |           |
  // 	  |  u(i-1)   |___________|           |
  // 	  |___________|           |   u(i+1)  |
  // 	  |           |           |___________|
  // ...|-----0-----|-----0-----|-----0-----|...
  // 	  |    i-1    |     i     |    i+1    |
  // 	  |+         -|+         -|+         -|
  //   i-3/2       i-1/2       i+1/2       i+3/2
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
  //                          -|
  //                          i+1/2
  //
  //
  //           |___________S0__________|
  //           |                       |
  //           |       |___________S1__________|
  //           |       |                       |
  //           |       |       |___________S2__________|
  //         ..|---o---|---o---|---o---|---o---|---o---|...
  //           | I{i-1}|  I{i} | I{i+1}| I{i+2}| I{i+3}|
  //                           |+
  //                          i+1/2
  //
  // WENO stencil: S{i} = [ I{i-2}, ..., I{i+3} ]
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  const auto epsilon  = static_cast<sc_t>(1e-6);
  const auto one	  = static_cast<sc_t>(1);
  const auto two	  = static_cast<sc_t>(2);
  const auto three	  = static_cast<sc_t>(3);
  const auto four	  = two*two;
  const auto five	  = three+two;
  const auto six	  = three*two;
  const auto seven	  = four+three;
  const auto ten      = five*two;
  const auto eleven   = five+six;
  const auto twelve   = six*two;
  const auto thirteen = six+seven;

  const auto oneOvfour= one/four;
  const auto oneOvsix = one/six;
  const auto oneOvten = one/ten;
  const auto threeOvten = three/ten;
  const auto sixOvten = six/ten;
  const auto thirteenOvtwelve = thirteen/twelve;

  {
    // u_{i+1/2}^{-}
    const auto p0 = (two*qim2 - seven*qim1 + eleven*qi)*oneOvsix;
    const auto p1 = ( -qim1   + five*qi    + two*qip1 )*oneOvsix;
    const auto p2 = (two*qi   + five*qip1  - qip2     )*oneOvsix;

    sc_t dp0_dq[5];
    sc_t dp1_dq[5];
    sc_t dp2_dq[5];
 
    dp0_dq[0] = 1./3.;
    dp0_dq[1] = -7./6.;
    dp0_dq[2] = 11./6.;
		dp0_dq[3] = 0.;
    dp0_dq[4] = 0.;

    dp1_dq[0] = 0.;
    dp1_dq[1] = -1./6.;
    dp1_dq[2] = 5./6.;
    dp1_dq[3] = 1./3.;
    dp1_dq[4] = 0.;

    dp2_dq[0] = 0.;
    dp2_dq[1] = 0.;
    dp2_dq[2] = 1./3.;
    dp2_dq[3] = 5./6.;
    dp2_dq[4] = -1./6.;



    const auto B0 = thirteenOvtwelve* std::pow(qim2-two*qim1+qi  , two) +
		    oneOvfour*std::pow(qim2-four*qim1+three*qi, two);

    const auto B1 = thirteenOvtwelve* std::pow(qim1 -two*qi +qip1, two) +
                    oneOvfour*std::pow(qim1-qip1, two);

    const auto B2 = thirteenOvtwelve* std::pow(qi  -two*qip1+qip2, two) +
		    oneOvfour*std::pow(three*qi-four*qip1+qip2, two);

    sc_t dB0_dq[5];
    sc_t dB1_dq[5];
    sc_t dB2_dq[5];

    dB0_dq[0] = (13.*(qim2-2.*qim1+qi))/6.+(qim2-4.*qim1+3.*qi)/2.;
    dB0_dq[1] = -(13.*(qim2-2.*qim1+qi))/3.-2.*(qim2-4.*qim1+3.*qi);
    dB0_dq[2] = (13.*(qim2-2.*qim1+qi))/6.+(3.*(qim2-4.*qim1+3.*qi))/2.;
    dB0_dq[3] = 0.;
    dB0_dq[4] = 0.;

    dB1_dq[0] = 0.;
    dB1_dq[1] = (13.*(qip1+qim1-2.*qi))/6.+(qim1-qip1)/2.;
    dB1_dq[2] = -(13.*(qip1+qim1-2.*qi))/3.;
    dB1_dq[3] = (13.*(qip1+qim1-2.*qi))/6.-(qim1-qip1)/2.;
    dB1_dq[4] = 0.;

    dB2_dq[0] = 0.;
    dB2_dq[1] = 0.;
    dB2_dq[2] = (13.*(qip2-2.*qip1+qi))/6.+(3.*(qip2-4.*qip1+3.*qi))/2.;
    dB2_dq[3] = -(13.*(qip2-2.*qip1+qi))/3.-2.*(qip2-4.*qip1+3.*qi);
    dB2_dq[4] = (13.*(qip2-2.*qip1+qi))/6.+(qip2-4.*qip1+3.*qi)/2.;
 

    const auto alpha0    = oneOvten/std::pow(epsilon   + B0, two);
    const auto alpha1    = sixOvten/std::pow(epsilon   + B1, two);
    const auto alpha2    = threeOvten/std::pow(epsilon + B2, two);
    const auto alphaSInv = one/(alpha0 + alpha1 + alpha2);

    sc_t dalpha0_dq[5]; // = dalpha0 /dB0 dB0/dq
    sc_t dalpha1_dq[5];
    sc_t dalpha2_dq[5];
    sc_t dalphaSInv_dq[5]; // = dalphaSinv/dalph0 dalph0/dq + ...

    for (int i = 0; i < 5; i++){
      dalpha0_dq[i] = -1./(5. * std::pow(B0 + epsilon,3.))*dB0_dq[i];
      dalpha1_dq[i] = -6./(5. * std::pow(B1 + epsilon,3.))*dB1_dq[i];
      dalpha2_dq[i] = -3./(5. * std::pow(B2 + epsilon,3.))*dB2_dq[i];
      dalphaSInv_dq[i] = -1./std::pow(alpha2+alpha1+alpha0,2.)*(dalpha0_dq[i] + dalpha1_dq[i] + dalpha2_dq[i]); 
    }

    const auto w0 = alpha0*alphaSInv;
    const auto w1 = alpha1*alphaSInv;
    const auto w2 = alpha2*alphaSInv;

    for (int i = 0; i < 5; i++){
      duNeg_dq[i] =  (dalpha0_dq[i]*alphaSInv + dalphaSInv_dq[i]*alpha0)*p0 + dp0_dq[i]*w0;
      duNeg_dq[i] += (dalpha1_dq[i]*alphaSInv + dalphaSInv_dq[i]*alpha1)*p1 + dp1_dq[i]*w1;
      duNeg_dq[i] += (dalpha2_dq[i]*alphaSInv + dalphaSInv_dq[i]*alpha2)*p2 + dp2_dq[i]*w2;
    }
    duNeg_dq[5] = 0.;
    uNeg = w0*p0 + w1*p1 + w2*p2;

  }

  {
    // u_{i+1/2}^{+}
    const auto p0 = (-qim1        + five*qi    + two*qip1 )*oneOvsix;
    const auto p1 = (two*qi       + five*qip1  - qip2     )*oneOvsix;
    const auto p2 = (eleven*qip1  - seven*qip2 + two*qip3 )*oneOvsix;

    sc_t dp0_dq[5];
    sc_t dp1_dq[5];
    sc_t dp2_dq[5];
 
    dp0_dq[0] = -1./6.;
    dp0_dq[1] = 5./6.;
    dp0_dq[2] = 1./3.;
		dp0_dq[3] = 0.;
    dp0_dq[4] = 0.;

    dp1_dq[0] = 0.;
    dp1_dq[1] = 1./3.;
    dp1_dq[2] = 5./6.;
    dp1_dq[3] = -1./6.;
    dp1_dq[4] = 0.;

    dp2_dq[0] = 0.;
    dp2_dq[1] = 0.;
    dp2_dq[2] = 11./6.;
    dp2_dq[3] = -7./6.;
    dp2_dq[4] = 1./3.;

    const auto B0 = thirteenOvtwelve* std::pow(qim1-two*qi+qip1   , two) +
		    oneOvfour*std::pow(qim1-four*qi+three*qip1, two);

    const auto B1 = thirteenOvtwelve* std::pow(qi  -two*qip1 +qip2, two) +
		    oneOvfour*std::pow(qi-qip2, two);

    const auto B2 = thirteenOvtwelve* std::pow(qip1-two*qip2+qip3 , two) +
		    oneOvfour*std::pow(three*qip1-four*qip2+qip3, two);

    sc_t dB0_dq[5];
    sc_t dB1_dq[5];
    sc_t dB2_dq[5];

    dB0_dq[0] = (3.*qip1+qim1-4.*qi)/2.+(13.*(qip1+qim1-2.*qi))/6.;
    dB0_dq[1] = -2.*(3.*qip1+qim1-4.*qi)-(13.*(qip1+qim1-2*qi))/3.;
    dB0_dq[2] = (3.*(3.*qip1+qim1-4.*qi))/2.+(13.*(qip1+qim1-2.*qi))/6.;
    dB0_dq[3] = 0.;
    dB0_dq[4] = 0.;

    dB1_dq[0] = 0.;
    dB1_dq[1] = (13.*(qip2-2*qip1+qi))/6.+(qi-qip2)/2.;
    dB1_dq[2] = -(13.*(qip2-2.*qip1+qi))/3.;
    dB1_dq[3] = (13.*(qip2-2.*qip1+qi))/6.-(qi-qip2)/2.;
    dB1_dq[4] = 0.;

    dB2_dq[0] = 0.;
    dB2_dq[1] = 0.;
    dB2_dq[2] = (13.*(qip3-2.*qip2+qip1))/6.+(3.*(qip3-4.*qip2+3.*qip1))/2.;
    dB2_dq[3] = -(13.*(qip3-2.*qip2+qip1))/3.-2.*(qip3-4.*qip2+3.*qip1);
    dB2_dq[4] = (13.*(qip3-2.*qip2+qip1))/6.+(qip3-4.*qip2+3.*qip1)/2.;


    const auto alpha0    = threeOvten/std::pow(epsilon + B0, two);
    const auto alpha1    = sixOvten/std::pow(epsilon   + B1, two);
    const auto alpha2    = oneOvten/std::pow(epsilon   + B2, two);
    const auto alphaSInv = one/(alpha0 + alpha1 + alpha2);

    sc_t dalpha0_dq[5]; // = dalpha0 /dB0 dB0/dq
    sc_t dalpha1_dq[5];
    sc_t dalpha2_dq[5];
    sc_t dalphaSInv_dq[5]; // = dalphaSinv/dalph0 dalph0/dq + ...

    for (int i = 0; i < 5; i++){
      dalpha0_dq[i] = -3./(5. * std::pow(B0 + epsilon,3.))*dB0_dq[i];
      dalpha1_dq[i] = -6./(5. * std::pow(B1 + epsilon,3.))*dB1_dq[i];
      dalpha2_dq[i] = -1./(5. * std::pow(B2 + epsilon,3.))*dB2_dq[i];
      dalphaSInv_dq[i] = -1./std::pow(alpha2+alpha1+alpha0,2.)*(dalpha0_dq[i] + dalpha1_dq[i] + dalpha2_dq[i]); 
    }

    const auto w0 = alpha0*alphaSInv;
    const auto w1 = alpha1*alphaSInv;
    const auto w2 = alpha2*alphaSInv;

    for (int i = 0; i < 5; i++){
      duPos_dq[i+1] =  (dalpha0_dq[i]*alphaSInv + dalphaSInv_dq[i]*alpha0)*p0 + dp0_dq[i]*w0;
      duPos_dq[i+1] += (dalpha1_dq[i]*alphaSInv + dalphaSInv_dq[i]*alpha1)*p1 + dp1_dq[i]*w1;
      duPos_dq[i+1] += (dalpha2_dq[i]*alphaSInv + dalphaSInv_dq[i]*alpha2)*p2 + dp2_dq[i]*w2;
    }
    duPos_dq[0] = 0.;
    uPos = w0*p0 + w1*p1 + w2*p2;

  }
}


}}
#endif
