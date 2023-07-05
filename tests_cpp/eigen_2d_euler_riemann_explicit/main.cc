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

  const auto probId  = pda::Euler2d::Riemann;
  auto appObj      = pda::create_problem_eigen(meshObj, probId, order, 2);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;

  state_t state(appObj.initialCondition());

  FomObserver<state_t> Obs("riemann2d_solution.bin", 100);
  const auto dt = 0.01;
  const auto Nsteps = pressio::ode::StepCount(0.6/dt);

  auto stepperObj = pressio::ode::create_ssprk3_stepper(appObj);
  pressio::ode::advance_n_steps(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
