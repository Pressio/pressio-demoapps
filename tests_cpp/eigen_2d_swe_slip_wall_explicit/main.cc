
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/swe2d.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

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
  auto appObj = pda::createProblemEigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;

  state_t state = appObj.initialCondition();

  const auto dt = 0.005;
  const auto Nsteps = 7.5/dt;
  FomObserver<state_t> Obs("swe_slipWall2d_solution.bin", 1);
  auto stepperObj = pressio::ode::create_rk4_stepper(state, appObj);
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0.0, dt, Nsteps, Obs);

  pressio::log::finalize();
  return 0;
}
