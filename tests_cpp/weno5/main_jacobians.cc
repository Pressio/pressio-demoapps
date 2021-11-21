#include <iostream>
#include <stdio.h>
#include <cmath>
#include <complex>
#include "pressiodemoapps/weno.hpp"
int main(){

  auto passedString = "PASS";
  double tol = 1e-12;
  double eps = 1e-18;
  int nVars = 4;
  using scalar_t = double;
  using complex_scalar_t = std::complex<double>;

  scalar_t U[6];
	U[0] = 0.4;
  U[1] = 8.;
  U[2] = 3.;
  U[3] = 0.5;
  U[4] = 0.9;
  U[5] = 128.3;

  scalar_t UNegGrad[6];
  scalar_t UPosGrad[6];
  scalar_t uNeg;
  scalar_t uPos;
  scalar_t uNeg_true;
  scalar_t uPos_true;

  ::pressiodemoapps::impl::weno5WithGrad(uNeg,uPos,UNegGrad,UPosGrad,U[0],U[1],U[2],U[3],U[4],U[5]);
  ::pressiodemoapps::impl::weno5(uNeg_true,uPos_true,U[0],U[1],U[2],U[3],U[4],U[5]);
  if (std::abs( uNeg_true - uNeg ) > tol){
       passedString = "FAILED";}
  if (std::abs( uPos_true - uPos ) > tol){
       passedString = "FAILED";}

  complex_scalar_t UCS[6];
	UCS[0] = 0.4;
  UCS[1] = 8.;
  UCS[2] = 3.;
  UCS[3] = 0.5;
  UCS[4] = 0.9;
  UCS[5] = 128.3;

  scalar_t UNegGradCS[6];
  scalar_t UPosGradCS[6];
  complex_scalar_t oneJ(0.,eps);
  for (int i = 0; i < 6; i++){
    complex_scalar_t uNegCS(0.,0.);
    complex_scalar_t uPosCS(0.,0.);
    UCS[i] += oneJ;
    ::pressiodemoapps::impl::weno5(uNegCS,uPosCS,UCS[0],UCS[1],UCS[2],UCS[3],UCS[4],UCS[5]);
    UNegGradCS[i] = 1./eps*std::imag(uNegCS);
    UPosGradCS[i] = 1./eps*std::imag(uPosCS);
    UCS[i] -= oneJ;

    if (std::abs( UNegGradCS[i]  - UNegGrad[i] ) > tol){
       passedString = "FAILED";
       std::cout << "Failed on uNeg " << i << std::endl;}

    if (std::abs( UPosGradCS[i]  - UPosGrad[i] ) > tol){
       passedString = "FAILED";
       std::cout << "Failed on uPos " << i << std::endl;}

  }
  std::cout << passedString << std::endl;
} 
