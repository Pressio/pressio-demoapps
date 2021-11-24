
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "../observer.hpp"
#include <random>

template<class scalar_type, class prim_t>
void initcond(int i, const scalar_type x, const scalar_type y, prim_t & prim)
{
  prim[0] = 1. + 0.1*x + 0.2*y;
  prim[1] = 0.1;
  prim[2] = 0.2;
  prim[3] = 0.5;
}

template<class T>
void writeToFileRank1(const T & obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (size_t i=0; i<obj.rows(); i++){
    file << std::setprecision(14) << obj(i) << " \n";
  }
  file.close();
}

template<class T>
void writeToFileSparseMat(const T & obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (size_t i=0; i<obj.rows(); i++){
    for (size_t j=0; j<obj.cols(); j++){
      file << std::setprecision(14) << obj.coeff(i,j) << " ";
    }
    file << " \n";
  }
  file.close();
}

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;

  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
#ifdef USE_WENO5
  const auto order   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#else
  const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probId  = pda::Euler2d::testingonlyneumann;
  auto appObj      = pda::create_problem_eigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;
  using scalar_t = typename app_t::scalar_type;

  const auto gamma = appObj.gamma();
  const auto stateSize = appObj.totalDofStencilMesh();
  state_t state(stateSize);
  std::array<scalar_t, 4> prim;
  const auto &x = meshObj.viewX();
  const auto &y = meshObj.viewY();
  for (int i=0; i<x.size(); ++i){
    initcond(i, x(i), y(i), prim);
    const auto ind = i*4;
    state(ind)   = prim[0];
    state(ind+1) = prim[0]*prim[1];
    state(ind+2) = prim[0]*prim[2];
    state(ind+3) = pressiodemoapps::eulerEquationsComputeEnergyFromPrimitive(gamma, prim);
  }

  auto time = 0.0;
  auto velo = appObj.createVelocity();
  auto jac  = appObj.createJacobian();
  appObj.velocityAndJacobian(state, time, velo, jac);
  //appObj.jacobian(state, time, jac);

  writeToFileRank1(state,   "state.txt");
  writeToFileRank1(velo,    "velo.txt");
  writeToFileSparseMat(jac, "jacobian.txt");

  return 0;
}
