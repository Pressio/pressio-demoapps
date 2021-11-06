
#include <iostream>
#include <stdio.h>
#include <cmath>
#include "pressiodemoapps/euler1d.hpp"

int main(){

  using scalar_t = double;
  auto passedString = "PASS";
  scalar_t eps = 1e-6;
  scalar_t tol = 1e-4;

  scalar_t gamma = 1.4;
  scalar_t rhoL = 1.;
  scalar_t uL = -1.;
  scalar_t pL = 1.;
  Eigen::Matrix<scalar_t,1,-1> UL(3);
  UL(0) = rhoL;
  UL(1) = rhoL*uL;
  UL(2) = pL / (gamma - 1.) + 0.5*rhoL*uL*uL;

  scalar_t rhoR = 0.125;
  scalar_t uR = -0.8;
  scalar_t pR = 0.1;
  Eigen::Matrix<scalar_t,1,-1> UR(3);
  UR(0) = rhoR;
  UR(1) = rhoR*uR;
  UR(2) = pR / (gamma - 1.) + 0.5*rhoR*uR*uR;


  Eigen::Matrix<scalar_t,1,-1> flux(3);
  Eigen::Matrix<scalar_t,1,-1> fluxBase(3);

  Eigen::Matrix<scalar_t,-1,-1> JL(3,3);
  Eigen::Matrix<scalar_t,-1,-1> JR(3,3);
  Eigen::Matrix<scalar_t,-1,-1> JL_FD(3,3);
  Eigen::Matrix<scalar_t,-1,-1> JR_FD(3,3);

  pressiodemoapps::ee::impl::eeRusanovFluxJacobianThreeDof(JL,JR,UL,UR,gamma);
  pressiodemoapps::ee::impl::eeRusanovFluxThreeDof(fluxBase,UL,UR,gamma);
  for (int i = 0;i < 3; i++){

    UL(i) += eps;
    pressiodemoapps::ee::impl::eeRusanovFluxThreeDof(flux,UL,UR,gamma);
    for (int j = 0 ; j < 3;j++){
      JL_FD(j,i) = 1./eps*(flux(j) - fluxBase(j) );
      if (std::abs( JL_FD(j,i) - JL(j,i) ) > tol){
        passedString = "FAILED";
        std::cout << "FAILED here " << i << " " << j << std::endl;}
    }
    UL(i) -= eps;

    UR(i) += eps;
    pressiodemoapps::ee::impl::eeRusanovFluxThreeDof(flux,UL,UR,gamma);
    for (int j = 0 ; j < 3;j++){
      JR_FD(j,i) = 1./eps*(flux(j) - fluxBase(j) );
      if (std::abs( JR_FD(j,i) - JR(j,i) ) > tol){
        passedString = "FAILED";}
    }
    UR(i) -= eps;

  }

  std::cout << passedString << std::endl;
}
