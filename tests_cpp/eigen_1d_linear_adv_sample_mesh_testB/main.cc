
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/advection.hpp"
#include "../observer.hpp"

template<class T>
void writeToFile(const T& obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (size_t i=0; i<obj.size(); i++){
    file << std::setprecision(14) << obj(i) << " \n";
  }
  file.close();
}

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
#ifdef USE_WENO5
  const auto order   = pda::InviscidFluxReconstruction::Weno5;
#else
  const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probid = pda::Advection1d::PeriodicLinear;
  auto appObj      = pda::createProblemEigen(meshObj, probid, order);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;
  using scalar_t = typename app_t::scalar_type;

  auto ic = appObj.initialCondition();
  state_t state(ic);

  auto time = 0.0;
  auto velo = appObj.createVelocity();
  appObj.velocity(state, time, velo);
  writeToFile(velo,  "velo.txt");
  writeToFile(state, "state.txt");

  // for (int i=0; i<velo.size(); ++i){
  //   std::cout << i << " " << i % 4 << " " << std::setprecision(14) << velo(i) << "\n";
  // }

  return 0;
}