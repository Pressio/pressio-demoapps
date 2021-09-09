
#include "pressio_ode_explicit.hpp"
#include "swe2d.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
{

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
#ifdef USE_WENO5
  const auto order   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#else
  const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probId  = pda::Swe2d::SlipWall;
  auto appObj      = pda::createProblemEigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using app_state_t = typename app_t::state_type;
  using ode_state_t = pressio::containers::Vector<app_state_t>;

  auto ic = appObj.initialCondition();
  ode_state_t state(ic);

  const auto dt = 0.005;
  const auto Nsteps = 7.5/dt;
  auto time = 0.0;
  FomObserver<ode_state_t> Obs("swe_slipWall2d_solution.bin", 1);


  auto stepperObj = pressio::ode::createRungeKutta4Stepper(state, appObj);

  pressio::ode::advanceNSteps(stepperObj, state, time, dt, Nsteps, Obs);
  return 0;
}
