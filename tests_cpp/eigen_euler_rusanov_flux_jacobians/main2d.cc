
#include <iostream>
#include <stdio.h>
#include <cmath>
#include "pressiodemoapps/euler2d.hpp"

int main(){
  auto passedString = "PASS";
  double tol = 1e-4;
  double eps = 1e-6;
  using scalar_t = double;
  scalar_t gamma = 1.4;
  scalar_t rhoL = 1.;
  scalar_t uL = -1.;
  scalar_t vL = 2.;
  scalar_t pL = 1.;

  Eigen::Matrix<scalar_t,1,-1> UL(4);
  UL(0) = rhoL;
  UL(1) = rhoL*uL;
  UL(2) = rhoL*vL;
  UL(3) = pL / (gamma - 1.) + 0.5*rhoL*uL*uL + 0.5*rhoL*vL*vL;

  scalar_t rhoR = 0.125;
  scalar_t uR = -0.8;
  scalar_t vR = 0.3;
  scalar_t pR = 0.1;
  Eigen::Matrix<scalar_t,1,-1> UR(4);
  UR(0) = rhoR;
  UR(1) = rhoR*uR;
  UR(2) = rhoR*vR;
  UR(3) = pR / (gamma - 1.) + 0.5*rhoR*uR*uR + 0.5*rhoR*vR*vR;
  scalar_t normals[2];
  normals[0] = 1;
  normals[1] = 0;
  Eigen::Matrix<scalar_t,1,-1> flux(4);
  Eigen::Matrix<scalar_t,1,-1> fluxBase(4);
  Eigen::Matrix<scalar_t,-1,-1> JL(4,4);
  Eigen::Matrix<scalar_t,-1,-1> JR(4,4);
  Eigen::Matrix<scalar_t,-1,-1> JL_FD(4,4);
  Eigen::Matrix<scalar_t,-1,-1> JR_FD(4,4);
  pressiodemoapps::ee::impl::ee_rusanov_flux_jacobian_four_dof(JL,JR,UL,UR,normals,gamma);
  pressiodemoapps::ee::impl::ee_rusanov_flux_four_dof(fluxBase,UL,UR,normals,gamma);
  for (int i = 0;i < 4; i++){
    UL(i) += eps;
    pressiodemoapps::ee::impl::ee_rusanov_flux_four_dof(flux,UL,UR,normals,gamma);
    for (int j = 0 ; j < 4;j++){
      JL_FD(j,i) = 1./eps*(flux(j) - fluxBase(j) );
      if (std::abs( JL_FD(j,i) - JL(j,i) ) > tol) passedString = "FAILED";
    }
    UL(i) -= eps;

    UR(i) += eps;
    pressiodemoapps::ee::impl::ee_rusanov_flux_four_dof(flux,UL,UR,normals,gamma);
    for (int j = 0 ; j < 4;j++){
      JR_FD(j,i) = 1./eps*(flux(j) - fluxBase(j) );
      if (std::abs( JR_FD(j,i) - JR(j,i) ) > tol) passedString = "FAILED";
    }
    UR(i) -= eps;
  }

  normals[0] = 0;
  normals[1] = 1;
  pressiodemoapps::ee::impl::ee_rusanov_flux_jacobian_four_dof(JL,JR,UL,UR,normals,gamma);
  pressiodemoapps::ee::impl::ee_rusanov_flux_four_dof(fluxBase,UL,UR,normals,gamma);
  for (int i = 0;i < 4; i++){
    UL(i) += eps;
    pressiodemoapps::ee::impl::ee_rusanov_flux_four_dof(flux,UL,UR,normals,gamma);
    for (int j = 0 ; j < 4;j++){
      JL_FD(j,i) = 1./eps*(flux(j) - fluxBase(j) );
      if (std::abs( JL_FD(j,i) - JL(j,i) ) > tol) passedString = "FAILED";
    }
    UL(i) -= eps;

    UR(i) += eps;
    pressiodemoapps::ee::impl::ee_rusanov_flux_four_dof(flux,UL,UR,normals,gamma);
    for (int j = 0 ; j < 4;j++){
      JR_FD(j,i) = 1./eps*(flux(j) - fluxBase(j) );
      if (std::abs( JR_FD(j,i) - JR(j,i) ) > tol) passedString = "FAILED";
    }
    UR(i) -= eps;
  }


  std::cout << passedString << std::endl;

}
