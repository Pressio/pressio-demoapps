
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/advection.hpp"
#include "../observer.hpp"

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
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
#ifdef USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#elif USE_WENO5
  const auto order   = pda::InviscidFluxReconstruction::Weno5;
#else
  const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probid = pda::Advection1d::PeriodicLinear;
  auto appObj      = pda::createProblemEigen(meshObj, probid, order);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;

  auto ic = appObj.initialCondition();
  state_t state(ic);

  auto time = 0.0;
  auto velo = appObj.createVelocity();
  auto jac  = appObj.createJacobian();
  appObj.velocity(state, time, velo);
  appObj.jacobian(state, time, jac);

  writeToFileRank1(state, "state.txt");
  writeToFileRank1(velo,  "velo.txt");
  writeToFileSparseMat(jac,   "jacobian.txt");

  return 0;
}
