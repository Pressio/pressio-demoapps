
#ifndef PRESSIODEMOAPPS_EE_JACOBIANS_HPP_
#define PRESSIODEMOAPPS_EE_JACOBIANS_HPP_

namespace pressiodemoapps{ namespace ee{ namespace impl{

template<class T, typename T2, typename sc_t>
void eeRusanovFluxJacobianThreeDof(T & JL, T & JR,
				   const T2 & qL,
				   const T2 & qR,
				   const sc_t gamma)
{

  constexpr sc_t one  = static_cast<sc_t>(1);
  constexpr sc_t two  = static_cast<sc_t>(2);
  constexpr sc_t three  = static_cast<sc_t>(3);
  constexpr sc_t half = one/two;
  const sc_t gm1 = gamma - one;
  constexpr sc_t es{1.e-30};

  // left
  const auto rL = qL(0);
  const auto uL = qL(1)/(rL + es);
  const auto pL = gm1*( qL(2) - half*rL*(uL*uL) );
  const auto HL = ( qL(2) + pL ) / rL;
  const auto aL = std::sqrt( gm1*(HL - half*(uL*uL)) );

  // right
  const auto rR = qR(0);
  const auto uR = qR(1)/(rR + es);
  const auto pR = gm1*( qR(2) - half*rR*(uR*uR) );
  const auto HR = ( qR(2) + pR ) / rR;
  const auto aR = std::sqrt( gm1*(HR - half*(uR*uR)) );

  // compute Roe average
  const auto r = std::sqrt(rR*rL);
  const auto RT = std::sqrt(rR/(rL));
  const auto u = (uL+RT*uR)/(one+RT);
  const auto H = ( HL+RT* HR)/(one+RT);
  const auto a = std::sqrt( gm1*(H - half*(u*u)) );
  const auto smax = std::abs(u) + a;


  // compute gradient of maximum eigen value
  sc_t gradL[3];
  sc_t gradR[3];

  const auto VMagSqrRoe = u*u + es;
  const auto VMagSqrL = uL*uL;
  const auto VVRoeL = uL*u ;
  const auto VMagSqrR = uR*uR;
  const auto VVRoeR = uR*u;
  const auto uRelRoe = u / std::sqrt(VMagSqrRoe);

  gradL[0] = one/(rL + r)*(-half*(uL + u)*uRelRoe  +  half*gm1/a * ( half*( VMagSqrRoe + VVRoeL) + half*(HL - H) - aL*aL / gm1 + half*(gamma - two)*VMagSqrL) );
  gradL[1] = one/(rL + r)*(uRelRoe - half*(gm1 * (u + gm1*uL) )/ (a));
  gradL[2] = half/(rL + r)*gamma*gm1 / (a);

  gradR[0] = one/(rR + r)*(-half*(uR + u)*uRelRoe +  half*gm1/a * ( half*( VMagSqrRoe + VVRoeR) + half*(HR - H) - aR*aR / gm1 + half*(gamma - two)*VMagSqrR) );
  gradR[1] = one/(rR + r)*(uRelRoe - half*(gm1 * (u + gm1*uR) )/ (a));
  gradR[2] = half/(rR + r)*gamma*gm1 / (a);

  // fill in Jacobians
  JL(0,0) = 0;
  JL(0,1) = half;
  JL(0,2) = 0;
  JL(1,0) = half*( half*gm1*uL*uL - uL*uL );
  JL(1,1) = half*( (three - gamma)*uL );
  JL(1,2) = half*gm1;
  JL(2,0) = half*( ( half*gm1*uL*uL - HL)*uL );
  JL(2,1) = half*( HL -  gm1*uL*uL );
  JL(2,2) = half*gamma*uL;

  for (int i = 0; i < 3; i++){
    JL(i,i) += half*smax;
    for (int j = 0; j < 3 ; j++){
      JL(i,j) += half*gradL[j]*(qL(i) - qR(i));
    }
  }

  JR(0,0) = 0;
  JR(0,1) = half;
  JR(0,2) = 0;
  JR(1,0) = half*( half*gm1*uR*uR - uR*uR );
  JR(1,1) = half*( (three - gamma)*uR );
  JR(1,2) = half*gm1;
  JR(2,0) = half*( ( half*gm1*uR*uR - HR)*uR );
  JR(2,1) = half*( HR -  gm1*uR*uR );
  JR(2,2) = half*gamma*uR;

  for (int i = 0; i < 3; i++){
    JR(i,i) -= half*smax;
    for (int j = 0; j < 3 ; j++){
      JR(i,j) += half*gradR[j]*(qL(i) - qR(i));
    }
  }
}


template<class T, typename T2, typename sc_t, typename normal_t>
void eeRusanovFluxJacobianFourDof(T & JL, T & JR,
         const T2 & qL,
         const T2 & qR,
        const normal_t & n,
         const sc_t gamma)
{
  constexpr auto one  = static_cast<sc_t>(1);
  constexpr auto two  = static_cast<sc_t>(2);
  constexpr auto half = one/two;
  const sc_t gm1 = gamma - one;
  constexpr sc_t es{1.e-30};

  // left
  const auto rL = qL(0);
  const auto uL = qL(1)/(rL + es);
  const auto vL = qL(2)/(rL + es);
  const auto unL = uL*n[0] + vL*n[1];
  const auto pL = gm1*( qL(3) - half*rL*(uL*uL+vL*vL) );
  const auto HL = ( qL(3) + pL ) / rL;
  const auto aL = std::sqrt( gm1*(HL - half*(uL*uL+vL*vL)) );

  // right
  const auto rR = qR(0);
  const auto uR = qR(1)/(rR + es);
  const auto vR = qR(2)/(rR + es);
  const auto unR = uR*n[0] + vR*n[1];
  const auto pR = gm1*( qR(3) - half*rR*(uR*uR+vR*vR) );
  const auto HR = ( qR(3) + pR ) / rR;
  const auto aR = std::sqrt( gm1*(HR - half*(uR*uR+vR*vR)) );

  // compute Roe average
  const auto r = std::sqrt(rR*rL);
  const auto RT = std::sqrt(rR/(rL));
  const auto u = (uL+RT*uR)/(one+RT);
  const auto v = (vL+RT*vR)/(one+RT);
  const auto H = ( HL+RT* HR)/(one+RT);
  const auto a = std::sqrt( gm1*(H - half*(u*u+v*v)) );
  const auto smax = std::sqrt(u*u+v*v) + a;

  // compute gradient of maximum eigen value
  sc_t gradL[4];
  sc_t gradR[4];

  const auto VMagSqrRoe = u*u + v*v + es;
  const auto VMagSqrL = uL*uL + vL*vL;
  const auto VVRoeL = uL*u + vL*v;
  const auto VMagSqrR = uR*uR + vR*vR;
  const auto VVRoeR = uR*u + vR*v;

  const auto uRelRoe = u / std::sqrt(VMagSqrRoe);
  const auto vRelRoe = v / std::sqrt(VMagSqrRoe);

  gradL[0] = one/(rL + r)*(-half*(uL + u)*uRelRoe - half*(vL + v)*vRelRoe +  half*gm1/a * ( half*( VMagSqrRoe + VVRoeL) + half*(HL - H) - aL*aL / gm1 + half*(gamma - two)*VMagSqrL) );
  gradL[1] = one/(rL + r)*(uRelRoe - half*(gm1 * (u + gm1*uL) )/ (a));
  gradL[2] = one/(rL + r)*(vRelRoe - half*(gm1 * (v + gm1*vL) )/ (a));
  gradL[3] = half/(rL + r)*gamma*gm1 / (a);

  gradR[0] = one/(rR + r)*(-half*(uR + u)*uRelRoe - half*(vR + v)*vRelRoe +  half*gm1/a * ( half*( VMagSqrRoe + VVRoeR) + half*(HR - H) - aR*aR / gm1 + half*(gamma - two)*VMagSqrR) );
  gradR[1] = one/(rR + r)*(uRelRoe - half*(gm1 * (u + gm1*uR) )/ (a));
  gradR[2] = one/(rR + r)*(vRelRoe - half*(gm1 * (v + gm1*vR) )/ (a));
  gradR[3] = half/(rR + r)*gamma*gm1 / (a);

  // fill in Jacobians
  JL(0,0) = 0;
  JL(0,1) = half*n[0];
  JL(0,2) = half*n[1];
  JL(0,3) = 0;

  JL(1,0) = half*( half*gm1*VMagSqrL*n[0] - uL*unL );
  JL(1,1) = half*( uL*n[0] - gm1*uL*n[0] + unL);
  JL(1,2) = half*( uL*n[1] - gm1*vL*n[0]);
  JL(1,3) = half*gm1*n[0];

  JL(2,0) = half*( half*gm1*VMagSqrL*n[1] - vL*unL );
  JL(2,1) = half*( vL*n[0] - gm1*uL*n[1]);
  JL(2,2) = half*( vL*n[1] - gm1*vL*n[1] + unL);
  JL(2,3) = half*gm1*n[1];

  JL(3,0) = half*( ( half*gm1*VMagSqrL - HL)*unL );
  JL(3,1) = half*( HL*n[0] -  gm1*uL*unL );
  JL(3,2) = half*( HL*n[1] - gm1*vL*unL );
  JL(3,3) = half*gamma*unL;

  for (int i = 0; i < 4; i++){
    JL(i,i) += half*smax;
    for (int j = 0; j < 4 ; j++){
      JL(i,j) += half*gradL[j]*(qL(i) - qR(i));
    }
  }

  JR(0,0) = 0;
  JR(0,1) = half*n[0];
  JR(0,2) = half*n[1];
  JR(0,3) = 0;

  JR(1,0) = half*( half*gm1*VMagSqrR*n[0] - uR*unR );
  JR(1,1) = half*( uR*n[0] - gm1*uR*n[0] + unR);
  JR(1,2) = half*( uR*n[1] - gm1*vR*n[0]);
  JR(1,3) = half*gm1*n[0];

  JR(2,0) = half*( half*gm1*VMagSqrR*n[1] - vR*unR );
  JR(2,1) = half*( vR*n[0] - gm1*uR*n[1]);
  JR(2,2) = half*( vR*n[1] - gm1*vR*n[1] + unR);
  JR(2,3) = half*gm1*n[1];

  JR(3,0) = half*( ( half*gm1*VMagSqrR - HR)*unR );
  JR(3,1) = half*( HR*n[0] -  gm1*uR*unR );
  JR(3,2) = half*( HR*n[1] - gm1*vR*unR );
  JR(3,3) = half*gamma*unR;

  for (int i = 0; i < 4; i++){
    JR(i,i) -= half*smax;
    for (int j = 0; j < 4 ; j++){
      JR(i,j) += half*gradR[j]*(qL(i) - qR(i));
    }
  }
}

template<class T, typename T2, typename sc_t, typename normal_t>
void eeRusanovFluxJacobianFiveDof(T & JL, T & JR,
				  const T2 & qL,
				  const T2 & qR,
				  const normal_t & n,
				  const sc_t gamma)
{

  constexpr auto one  = static_cast<sc_t>(1);
  constexpr auto two  = static_cast<sc_t>(2);
  constexpr auto half = one/two;
  const sc_t gm1 = gamma - one;
  constexpr sc_t es{1.e-30};

  // left
  const auto rL = qL(0);
  const auto uL = qL(1)/(rL + es);
  const auto vL = qL(2)/(rL + es);
  const auto wL = qL(3)/(rL + es);
  const auto unL = uL*n[0] + vL*n[1] + wL*n[2];
  const auto pL = gm1*( qL(4) - half*rL*(uL*uL+vL*vL+wL*wL) );
  const auto HL = ( qL(4) + pL ) / rL;
  const auto aL = std::sqrt( gm1*(HL - half*(uL*uL + vL*vL + wL*wL)) );

  // right
  const auto rR = qR(0);
  const auto uR = qR(1)/(rR + es);
  const auto vR = qR(2)/(rR + es);
  const auto wR = qR(3)/(rR + es);
  const auto unR = uR*n[0] + vR*n[1] + wR*n[2];
  const auto pR = gm1*( qR(4) - half*rR*(uR*uR+vR*vR+wR*wR) );
  const auto HR = ( qR(4) + pR ) / rR;
  const auto aR = std::sqrt( gm1*(HR - half*(uR*uR + vR*vR + wR*wR)) );

  // compute Roe average
  const auto r = std::sqrt(rR*rL);
  const auto RT = std::sqrt(rR/(rL));
  const auto u = (uL+RT*uR)/(one+RT);
  const auto v = (vL+RT*vR)/(one+RT);
  const auto w = (wL+RT*wR)/(one+RT);
  const auto H = ( HL+RT* HR)/(one+RT);
  const auto a = std::sqrt( gm1*(H - half*(u*u + v*v + w*w)) );
  const auto smax = std::sqrt(u*u+v*v+w*w) + a;
  // compute gradient of maximum eigen value
  sc_t gradL[5];
  sc_t gradR[5];

  const auto VMagSqrRoe = u*u + v*v + w*w + es;
  const auto VMagSqrL = uL*uL + vL*vL + wL*wL;
  const auto VVRoeL = uL*u + vL*v + wL*w;
  const auto VMagSqrR = uR*uR + vR*vR + wR*wR;
  const auto VVRoeR = uR*u + vR*v + wR*w;

  const auto uRelRoe = u / std::sqrt(VMagSqrRoe);
  const auto vRelRoe = v / std::sqrt(VMagSqrRoe);
  const auto wRelRoe = w / std::sqrt(VMagSqrRoe);

  gradL[0] = one/(rL + r)*(-half*(uL + u)*uRelRoe - half*(vL + v)*vRelRoe - half*(wL + w)*wRelRoe +  half*gm1/a * ( half*( VMagSqrRoe + VVRoeL) + half*(HL - H) - aL*aL / gm1 + half*(gamma - two)*VMagSqrL) );
  gradL[1] = one/(rL + r)*(uRelRoe - half*(gm1 * (u + gm1*uL) )/ (a));
  gradL[2] = one/(rL + r)*(vRelRoe - half*(gm1 * (v + gm1*vL) )/ (a));
  gradL[3] = one/(rL + r)*(wRelRoe - half*(gm1 * (w + gm1*wL) )/ (a));
  gradL[4] = half/(rL + r)*gamma*gm1 / (a);

  gradR[0] = one/(rR + r)*(-half*(uR + u)*uRelRoe - half*(vR + v)*vRelRoe - half*(wR + w)*wRelRoe +  half*gm1/a * ( half*( VMagSqrRoe + VVRoeR) + half*(HR - H) - aR*aR / gm1 + half*(gamma - two)*VMagSqrR) );
  gradR[1] = one/(rR + r)*(uRelRoe - half*(gm1 * (u + gm1*uR) )/ (a));
  gradR[2] = one/(rR + r)*(vRelRoe - half*(gm1 * (v + gm1*vR) )/ (a));
  gradR[3] = one/(rR + r)*(wRelRoe - half*(gm1 * (w + gm1*wR) )/ (a));
  gradR[4] = half/(rR + r)*gamma*gm1 / (a);

  JL(0,0) = 0;
  JL(0,1) = half*n[0];
  JL(0,2) = half*n[1];
  JL(0,3) = half*n[2];
  JL(0,4) = 0;

  JL(1,0) = half*( half*gm1*VMagSqrL*n[0] - uL*unL );
  JL(1,1) = half*( uL*n[0] - gm1*uL*n[0] + unL);
  JL(1,2) = half*( uL*n[1] - gm1*vL*n[0]);
  JL(1,3) = half*( uL*n[2] - gm1*wL*n[0]);
  JL(1,4) = half*gm1*n[0];

  JL(2,0) = half*( half*gm1*VMagSqrL*n[1] - vL*unL );
  JL(2,1) = half*( vL*n[0] - gm1*uL*n[1]);
  JL(2,2) = half*( vL*n[1] - gm1*vL*n[1] + unL);
  JL(2,3) = half*( vL*n[2] - gm1*wL*n[1]);
  JL(2,4) = half*gm1*n[1];

  JL(3,0) = half*( half*gm1*VMagSqrL*n[2] - wL*unL );
  JL(3,1) = half*( wL*n[0] - gm1*uL*n[2]);
  JL(3,2) = half*( wL*n[1] - gm1*vL*n[2]);
  JL(3,3) = half*( wL*n[2] - gm1*wL*n[2] + unL);
  JL(3,4) = half*gm1*n[2];

  JL(4,0) = half*( ( half*gm1*VMagSqrL - HL)*unL );
  JL(4,1) = half*( HL*n[0] -  gm1*uL*unL );
  JL(4,2) = half*( HL*n[1] -  gm1*vL*unL );
  JL(4,3) = half*( HL*n[2] -  gm1*wL*unL );
  JL(4,4) = half*gamma*unL;

  for (int i = 0; i < 5; i++){
    JL(i,i) += half*smax;
    for (int j = 0; j < 5 ; j++){
      JL(i,j) += half*gradL[j]*(qL(i) - qR(i));
    }
  }

  JR(0,0) = 0;
  JR(0,1) = half*n[0];
  JR(0,2) = half*n[1];
  JR(0,3) = half*n[2];
  JR(0,4) = 0;

  JR(1,0) = half*( half*gm1*VMagSqrR*n[0] - uR*unR );
  JR(1,1) = half*( uR*n[0] - gm1*uR*n[0] + unR);
  JR(1,2) = half*( uR*n[1] - gm1*vR*n[0]);
  JR(1,3) = half*( uR*n[2] - gm1*wR*n[0]);
  JR(1,4) = half*gm1*n[0];

  JR(2,0) = half*( half*gm1*VMagSqrR*n[1] - vR*unR );
  JR(2,1) = half*( vR*n[0] - gm1*uR*n[1]);
  JR(2,2) = half*( vR*n[1] - gm1*vR*n[1] + unR);
  JR(2,3) = half*( vR*n[2] - gm1*wR*n[1]);
  JR(2,4) = half*gm1*n[1];

  JR(3,0) = half*( half*gm1*VMagSqrR*n[2] - wR*unR );
  JR(3,1) = half*( wR*n[0] - gm1*uR*n[2]);
  JR(3,2) = half*( wR*n[1] - gm1*vR*n[2]);
  JR(3,3) = half*( wR*n[2] - gm1*wR*n[2] + unR);
  JR(3,4) = half*gm1*n[2];

  JR(4,0) = half*( ( half*gm1*VMagSqrR - HR)*unR );
  JR(4,1) = half*( HR*n[0] -  gm1*uR*unR );
  JR(4,2) = half*( HR*n[1] -  gm1*vR*unR );
  JR(4,3) = half*( HR*n[2] -  gm1*wR*unR );
  JR(4,4) = half*gamma*unR;

  for (int i = 0; i < 5; i++){
    JR(i,i) -= half*smax;
    for (int j = 0; j < 5 ; j++){
      JR(i,j) += half*gradR[j]*(qL(i) - qR(i));
    }
  }
}

}}}
#endif
