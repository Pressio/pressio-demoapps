
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "../observer.hpp"

int main()
{
  // pressio::log::initialize(pressio::logto::terminal);
  // pressio::log::setVerbosity({pressio::log::level::debug});

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");

#ifdef USE_WENO5
  constexpr auto order   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  constexpr auto order   = pda::InviscidFluxReconstruction::Weno3;
#else
  constexpr auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probId  = pda::Euler2d::PeriodicSmooth;
  auto appObj = pda::create_problem_eigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using state_t	= typename app_t::state_type;
  using jacob_t	= typename app_t::jacobian_type;

  state_t state = appObj.initialCondition();

  auto stepperObj = pressio::ode::create_implicit_stepper(
    pressio::ode::StepScheme::BDF1, appObj);

  using lin_solver_t = pressio::linearsolvers::Solver<
    pressio::linearsolvers::iterative::Bicgstab, jacob_t>;
  lin_solver_t linSolverObj;
  auto NonLinSolver=pressio::create_newton_solver(stepperObj, linSolverObj);
  NonLinSolver.setStopTolerance(1e-10);

  FomObserver<state_t> Obs("smooth2d_solution.bin", 5);

  const auto dt = 0.05;
  const auto Nsteps = pressio::ode::StepCount(40);
  pressio::ode::advance_n_steps(stepperObj, state, 0., dt, Nsteps, Obs, NonLinSolver);

  return 0;
}
