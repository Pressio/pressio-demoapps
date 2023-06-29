
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/advection_diffusion_reaction2d.hpp"
#include <random>

std::mt19937 gen(2195884);
std::uniform_real_distribution<double> distrib(1.1, 2.0);

template<class T>
void writeToFileRank1(const T & obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (int i=0; i<obj.rows(); i++){
    file << std::setprecision(14) << obj(i) << " \n";
  }
  file.close();
}

template<class T>
void writeToFileSparseMat(const T & obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (int i=0; i<obj.rows(); i++){
    for (int j=0; j<obj.cols(); j++){
      file << std::setprecision(14) << obj.coeff(i,j) << " ";
    }
    file << " \n";
  }
  file.close();
}

int main()
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
#ifdef USE_WENO5
  const auto inviscidScheme = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  const auto inviscidScheme= pda::InviscidFluxReconstruction::Weno3;
#else
  const auto inviscidScheme= pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probId  = pda::AdvectionDiffusionReaction2d::ProblemA;
  auto appObj = pda::create_problem_eigen(meshObj, probId, inviscidScheme);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;

  const auto stateSize = appObj.totalDofStencilMesh();
  state_t state(stateSize);
  const auto &x = meshObj.viewX();
  const auto &y = meshObj.viewY();
  for (int i=0; i<x.size(); ++i){
    const auto ind = i;
    state(ind) = 1. + 0.1*x(i) + 0.2*y(i);
  }

  auto time = 0.0;
  auto velo = appObj.createRightHandSide();
  auto jac  = appObj.createJacobian();

  appObj.rightHandSideAndJacobian(state, time, velo, jac);
  writeToFileRank1(state,   "state.txt");
  writeToFileRank1(velo,    "velo.txt");
  writeToFileSparseMat(jac, "jacobian.txt");
  return 0;
}
