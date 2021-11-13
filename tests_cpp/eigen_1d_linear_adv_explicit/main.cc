
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/advection.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
#ifdef USE_WENO5
  constexpr auto scheme   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  constexpr auto scheme   = pda::InviscidFluxReconstruction::Weno3;
#else
  constexpr auto scheme   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probid = pda::Advection1d::PeriodicLinear;
  auto appObj      = pda::createExplicitProblemEigen(meshObj, probid, scheme);

  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;
  state_t state = appObj.initialCondition();

  auto stepperObj = pressio::ode::create_rk4_stepper(state, appObj);
  FomObserver<state_t> Obs("sod1d_solution.bin", 10);

  const auto dt = 0.001;
  const auto Nsteps = 100;
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
