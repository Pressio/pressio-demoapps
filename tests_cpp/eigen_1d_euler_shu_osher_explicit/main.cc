
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler1d.hpp"
#include "../observer.hpp"

int main()
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
#ifdef USE_WENO5
  constexpr auto order   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  constexpr auto order   = pda::InviscidFluxReconstruction::Weno3;
#else
  constexpr auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  auto appObj      = pda::create_problem_eigen(meshObj, pda::Euler1d::ShuOsher, order);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;
  state_t state = appObj.initialCondition();

  auto stepperObj = pressio::ode::create_ssprk3_stepper(appObj);
  FomObserver<state_t> Obs("shu_osher1d_solution.bin", 1);

  const auto dt = 0.001;
  const auto Nsteps = pressio::ode::StepCount(1.8/dt);
  pressio::ode::advance_n_steps(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
