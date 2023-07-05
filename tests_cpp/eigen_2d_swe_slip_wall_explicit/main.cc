
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/swe2d.hpp"
#include "../observer.hpp"

int main()
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::info});

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
#ifdef USE_WENO5
  const auto order   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#else
  const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probId  = pda::Swe2d::SlipWall;
  auto appObj = pda::create_problem_eigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;

  state_t state = appObj.initialCondition();

  const auto dt = 0.01;
  const auto Nsteps = pressio::ode::StepCount(2./dt);
  FomObserver<state_t> Obs("swe_slipWall2d_solution.bin", 10);
  auto stepperObj = pressio::ode::create_rk4_stepper(appObj);
  pressio::ode::advance_n_steps(stepperObj, state, 0.0, dt, Nsteps, Obs);

  pressio::log::finalize();
  return 0;
}
