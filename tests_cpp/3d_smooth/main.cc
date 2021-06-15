
#include "pressio_ode_explicit.hpp"
#include "euler3d.hpp"
#include "../observer.hpp"

template<class T>
void writeToFile(const T& obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (size_t i=0; i<obj.extent(0); i++){
    file << std::setprecision(14) << obj(i) << " \n";
  }
  file.close();
}

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
  const auto order   = pda::ReconstructionType::firstOrder;

  const auto probId  = pda::Euler3d::PeriodicSmooth;
  auto appObj      = pda::createEuler3dEigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using app_state_t = typename app_t::state_type;
  using ode_state_t = pressio::containers::Vector<app_state_t>;
  using scalar_t = typename app_t::scalar_type;

  const auto stateSize = appObj.totalDofStencilMesh();
  ode_state_t state(appObj.initialCondition());
  // writeToFile(state,  "ic.txt");

  auto stepperObj = pressio::ode::createRungeKutta4Stepper(state, appObj);

  const auto dt = 0.001;
  const auto Nsteps = 2./dt;
  FomObserver<ode_state_t> Obs("solution.bin", 100);
  pressio::ode::advanceNSteps(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
