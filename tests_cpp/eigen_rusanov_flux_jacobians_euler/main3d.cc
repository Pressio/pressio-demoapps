
#include <iostream>
#include <stdio.h>
#include <cmath>
#include "pressiodemoapps/euler3d.hpp"

int main(){
  auto passedString = "PASS";
  double tol = 1e-4;
  double eps = 1e-6;
  using scalar_t = double;
  scalar_t gamma = 1.4;
  scalar_t rhoL = 1.;
  scalar_t uL = -1.;
  scalar_t vL = 2.;
  scalar_t wL = -0.3;
  scalar_t pL = 1.;
  Eigen::Matrix<scalar_t,1,-1> UL(5);
  UL(0) = rhoL;
  UL(1) = rhoL*uL;
  UL(2) = rhoL*vL;
  UL(3) = rhoL*wL;
  UL(4) = pL / (gamma - 1.) + 0.5*rhoL*uL*uL + 0.5*rhoL*vL*vL + 0.5*rhoL*wL*wL;

  scalar_t rhoR = 0.125;
  scalar_t uR = -0.8;
  scalar_t vR = 0.3;
  scalar_t wR = 2.3;
  scalar_t pR = 0.1;
  Eigen::Matrix<scalar_t,1,-1> UR(5);
  UR(0) = rhoR;
  UR(1) = rhoR*uR;
  UR(2) = rhoR*vR;
  UR(3) = rhoR*wR;
  UR(4) = pR / (gamma - 1.) + 0.5*rhoR*uR*uR + 0.5*rhoR*vR*vR + 0.5*rhoR*wR*wR;

  scalar_t normals[3];
  normals[0] = 1;
  normals[1] = 0;
  normals[2] = 0;

  Eigen::Matrix<scalar_t,1,-1> flux(5);
  Eigen::Matrix<scalar_t,1,-1> fluxBase(5);
  Eigen::Matrix<scalar_t,-1,-1> JL(5,5);
  Eigen::Matrix<scalar_t,-1,-1> JR(5,5);
  Eigen::Matrix<scalar_t,-1,-1> JL_FD(5,5);
  Eigen::Matrix<scalar_t,-1,-1> JR_FD(5,5);
  pressiodemoapps::ee::impl::ee_rusanov_flux_jacobian_five_dof(JL,JR,UL,UR,normals,gamma);
  pressiodemoapps::ee::impl::ee_rusanov_flux_five_dof(fluxBase,UL,UR,normals,gamma);
  for (int i = 0;i < 5; i++){
    UL(i) += eps;
    pressiodemoapps::ee::impl::ee_rusanov_flux_five_dof(flux,UL,UR,normals,gamma);
    for (int j = 0 ; j < 5;j++){
      JL_FD(j,i) = 1./eps*(flux(j) - fluxBase(j) );
      if (std::abs( JL_FD(j,i) - JL(j,i) ) > tol){
        passedString = "FAILED";}
    }
    UL(i) -= eps;

    UR(i) += eps;
    pressiodemoapps::ee::impl::ee_rusanov_flux_five_dof(flux,UL,UR,normals,gamma);
    for (int j = 0 ; j < 5;j++){
      JR_FD(j,i) = 1./eps*(flux(j) - fluxBase(j) );
      if (std::abs( JR_FD(j,i) - JR(j,i) ) > tol){
        passedString = "FAILED";}
    }
    UR(i) -= eps;
  }

  normals[0] = 0;
  normals[1] = 1;
  normals[2] = 0;

  pressiodemoapps::ee::impl::ee_rusanov_flux_jacobian_five_dof(JL,JR,UL,UR,normals,gamma);
  pressiodemoapps::ee::impl::ee_rusanov_flux_five_dof(fluxBase,UL,UR,normals,gamma);
  for (int i = 0;i < 5; i++){
    UL(i) += eps;
    pressiodemoapps::ee::impl::ee_rusanov_flux_five_dof(flux,UL,UR,normals,gamma);
    for (int j = 0 ; j < 5;j++){
      JL_FD(j,i) = 1./eps*(flux(j) - fluxBase(j) );
      if (std::abs( JL_FD(j,i) - JL(j,i) ) > tol) passedString = "FAILED";
    }
    UL(i) -= eps;

    UR(i) += eps;
    pressiodemoapps::ee::impl::ee_rusanov_flux_five_dof(flux,UL,UR,normals,gamma);
    for (int j = 0 ; j < 5;j++){
      JR_FD(j,i) = 1./eps*(flux(j) - fluxBase(j) );
      if (std::abs( JR_FD(j,i) - JR(j,i) ) > tol) passedString = "FAILED";
    }
    UR(i) -= eps;
  }

  normals[0] = 0;
  normals[1] = 0;
  normals[2] = 1;

  pressiodemoapps::ee::impl::ee_rusanov_flux_jacobian_five_dof(JL,JR,UL,UR,normals,gamma);
  pressiodemoapps::ee::impl::ee_rusanov_flux_five_dof(fluxBase,UL,UR,normals,gamma);
  for (int i = 0;i < 5; i++){
    UL(i) += eps;
    pressiodemoapps::ee::impl::ee_rusanov_flux_five_dof(flux,UL,UR,normals,gamma);
    for (int j = 0 ; j < 5;j++){
      JL_FD(j,i) = 1./eps*(flux(j) - fluxBase(j) );
      if (std::abs( JL_FD(j,i) - JL(j,i) ) > tol) passedString = "FAILED";
    }
    UL(i) -= eps;

    UR(i) += eps;
    pressiodemoapps::ee::impl::ee_rusanov_flux_five_dof(flux,UL,UR,normals,gamma);
    for (int j = 0 ; j < 5;j++){
      JR_FD(j,i) = 1./eps*(flux(j) - fluxBase(j) );
      if (std::abs( JR_FD(j,i) - JR(j,i) ) > tol) passedString = "FAILED";
    }
    UR(i) -= eps;
  }


  std::cout << passedString << std::endl;

}
