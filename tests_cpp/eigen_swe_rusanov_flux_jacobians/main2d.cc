#include <Eigen/Dense>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include "pressiodemoapps/impl/swe_flux_jacobian.hpp"
#include "pressiodemoapps/impl/swe_fluxes.hpp"
#include <complex>
int main(){
  auto passedString = "PASS";
  double tol = 1e-12;
  double eps = 1e-15;

  int nVars = 3;
  using scalar_t = std::complex<double>;

  scalar_t g = 9.8;
	scalar_t hL = 0.4;
  scalar_t uL = 8.;
  scalar_t vL = 3.;

  Eigen::Matrix<scalar_t,1,-1> UL(nVars);
  UL(0) = hL;
  UL(1) = hL*uL;
  UL(2) = hL*vL;

  scalar_t hR = 0.125;
  scalar_t uR = -0.8;
  scalar_t vR = -0.3;
  Eigen::Matrix<scalar_t,1,-1> UR(nVars);
  UR(0) = hR;
  UR(1) = hR*uR;
  UR(2) = hR*vR;
  scalar_t normals[2];
  normals[0] = 1;
  normals[1] = 0;
  Eigen::Matrix<scalar_t,1,-1> flux(nVars);
  Eigen::Matrix<scalar_t,1,-1> fluxBase(nVars);
  Eigen::Matrix<scalar_t,-1,-1> JL(nVars,nVars);
  Eigen::Matrix<scalar_t,-1,-1> JR(nVars,nVars);
  Eigen::Matrix<scalar_t,-1,-1> JL_FD(nVars,nVars);
  Eigen::Matrix<scalar_t,-1,-1> JR_FD(nVars,nVars);
  pressiodemoapps::implswe::sweRusanovFluxJacobianThreeDof(JL,JR,UL,UR,normals,g);
  pressiodemoapps::implswe::sweRusanovFluxThreeDof(fluxBase,UL,UR,normals,g);
  const scalar_t oneJ(0.0,eps); 
  for (int i = 0;i < nVars; i++){
    UL(i) += oneJ;
    pressiodemoapps::implswe::sweRusanovFluxThreeDof(flux,UL,UR,normals,g);
    for (int j = 0 ; j < nVars;j++){
      JL_FD(j,i) = 1./eps*(std::imag(flux(j)));
      if (std::abs( JL_FD(j,i) - JL(j,i) ) > tol) passedString = "FAILED";
    }
    UL(i) -= oneJ;

    UR(i) += oneJ;
    pressiodemoapps::implswe::sweRusanovFluxThreeDof(flux,UL,UR,normals,g);
    for (int j = 0 ; j < nVars;j++){
      JR_FD(j,i) = 1./eps*std::imag(flux(j) );
      if (std::abs( JR_FD(j,i) - JR(j,i) ) > tol) passedString = "FAILED";
    }
    UR(i) -= oneJ;
  }

  normals[0] = 0;
  normals[1] = 1;
  pressiodemoapps::implswe::sweRusanovFluxJacobianThreeDof(JL,JR,UL,UR,normals,g);
  pressiodemoapps::implswe::sweRusanovFluxThreeDof(fluxBase,UL,UR,normals,g);
  for (int i = 0;i < nVars; i++){
    UL(i) += oneJ;
    pressiodemoapps::implswe::sweRusanovFluxThreeDof(flux,UL,UR,normals,g);
    for (int j = 0 ; j < nVars;j++){
      JL_FD(j,i) = 1./eps*std::imag(flux(j) );
      if (std::abs( JL_FD(j,i) - JL(j,i) ) > tol) passedString = "FAILED";
    }
    UL(i) -= oneJ;

    UR(i) += oneJ;
    pressiodemoapps::implswe::sweRusanovFluxThreeDof(flux,UL,UR,normals,g);
    for (int j = 0 ; j < nVars;j++){
      JR_FD(j,i) = 1./eps*std::imag(flux(j));
      if (std::abs( JR_FD(j,i) - JR(j,i) ) > tol) passedString = "FAILED";
    }
    UR(i) -= oneJ;
  }


  std::cout << passedString << std::endl;

}
