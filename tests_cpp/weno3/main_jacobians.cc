#include <iostream>
#include <stdio.h>
#include <cmath>
#include <complex>
#include "pressiodemoapps/weno.hpp"

int main(){

  auto passedString = "PASS";
  double tol = 1e-12;
  double eps = 1e-18;

  using scalar_t = double;
  using complex_scalar_t = std::complex<double>;

  scalar_t U[4];
	U[0] = 0.4;
  U[1] = 8.;
  U[2] = 3.;
  U[3] = 0.5;

  scalar_t uNegGrad[4];
  scalar_t uPosGrad[4];
  scalar_t uNeg;
  scalar_t uNeg_true;
  scalar_t uPos;
  scalar_t uPos_true;

  ::pressiodemoapps::impl::weno3(uNeg_true,uPos_true,U[0],U[1],U[2],U[3]);
  ::pressiodemoapps::impl::weno3WithGrad(uNeg,uPos,uNegGrad,uPosGrad,U[0],U[1],U[2],U[3]);

  if (std::abs( uNeg_true - uNeg ) > tol){
       passedString = "FAILED";}
  if (std::abs( uPos_true - uPos ) > tol){
       passedString = "FAILED";}


  complex_scalar_t UCS[4];
	UCS[0] = 0.4;
  UCS[1] = 8.;
  UCS[2] = 3.;
  UCS[3] = 0.5;

  scalar_t uNegGradCS[4];
  scalar_t uPosGradCS[4];
  complex_scalar_t oneJ(0.,eps);
  for (int i = 0; i < 4; i++){
    complex_scalar_t uNegCS(0.,0.);
    complex_scalar_t uPosCS(0.,0.);
    UCS[i] += oneJ;
    ::pressiodemoapps::impl::weno3(uNegCS,uPosCS,UCS[0],UCS[1],UCS[2],UCS[3]);
    uNegGradCS[i] = 1./eps*std::imag(uNegCS);
    uPosGradCS[i] = 1./eps*std::imag(uPosCS);
    UCS[i] -= oneJ;

    if (std::abs( uNegGradCS[i]  - uNegGrad[i] ) > tol){
       passedString = "FAILED";
       std::cout << "Failed on uNeg " << i << std::endl;}

    if (std::abs( uPosGradCS[i]  - uPosGrad[i] ) > tol){
       passedString = "FAILED";
       std::cout << "Failed on uPos " << i << std::endl;}

  }
  std::cout << passedString << std::endl;
} 
