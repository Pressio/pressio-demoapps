
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler3d.hpp"
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
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
#ifdef USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#else
  const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probId  = pda::Euler3d::PeriodicSmooth;
  auto appObj      = pda::create_problem_eigen(meshObj, probId, order);
  appObj.disableJacobian();
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;

  state_t state(appObj.initialCondition());
  auto stepperObj = pressio::ode::create_rk4_stepper(state, appObj);

  const auto dt = 0.001;
  const auto Nsteps = 2./dt;
  FomObserver<state_t> Obs("solution.bin", 100);
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
