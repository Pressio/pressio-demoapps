
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/advection1d.hpp"
#include "../observer.hpp"

template<class T>
void writeToFileRank1(const T & obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (decltype(obj.rows()) i=0; i<obj.rows(); i++){
    file << std::setprecision(14) << obj(i) << " \n";
  }
  file.close();
}

template<class T>
void writeToFileSparseMat(const T & obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (decltype(obj.rows()) i=0; i<obj.rows(); i++){
    for (decltype(obj.cols()) j=0; j<obj.cols(); j++){
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
#ifdef USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#elif USE_WENO5
  const auto order   = pda::InviscidFluxReconstruction::Weno5;
#else
  const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probid = pda::Advection1d::PeriodicLinear;
  auto appObj      = pda::create_problem_eigen(meshObj, probid, order);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;

  auto ic = appObj.initialCondition();
  state_t state(ic);

  auto time = 0.0;
  auto velo = appObj.createRightHandSide();
  auto jac  = appObj.createJacobian();
  appObj.rightHandSide(state, time, velo);
  appObj.jacobian(state, time, jac);

  writeToFileRank1(state, "state.txt");
  writeToFileRank1(velo,  "velo.txt");
  writeToFileSparseMat(jac,   "jacobian.txt");

  return 0;
}
