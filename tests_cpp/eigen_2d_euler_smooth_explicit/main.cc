
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "../observer.hpp"

int main()
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

  const auto probId  = pda::Euler2d::PeriodicSmooth;
  auto appObj      = pda::create_problem_eigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;
  state_t state(appObj.initialCondition());

  auto stepperObj = pressio::ode::create_rk4_stepper(appObj);
  FomObserver<state_t> Obs("smooth2d_solution.bin", 5);

  const auto dt = 0.001;
  const auto Nsteps = pressio::ode::StepCount(2.0/dt);
  pressio::ode::advance_n_steps(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
