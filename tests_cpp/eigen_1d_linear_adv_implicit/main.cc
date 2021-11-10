
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/advection.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");

#ifdef USE_WENO5
  constexpr auto order   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  constexpr auto order   = pda::InviscidFluxReconstruction::Weno3;
#else
  constexpr auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probid = pda::Advection1d::PeriodicLinear;
  auto appObj      = pda::createProblemEigen(meshObj, probid, order);
  using app_t = decltype(appObj);
  using scalar_t	= typename app_t::scalar_type;
  using state_t	= typename app_t::state_type;
  using jacob_t	= typename app_t::jacobian_type;

  state_t state = appObj.initialCondition();

  auto stepperObj = pressio::ode::create_implicit_stepper(
    pressio::ode::StepScheme::BDF1, state, appObj);

  using lin_solver_t = pressio::linearsolvers::Solver<
    pressio::linearsolvers::iterative::Bicgstab, jacob_t>;
  lin_solver_t linSolverObj;

  auto NonLinSolver=
    pressio::nonlinearsolvers::create_newton_raphson(stepperObj, state, linSolverObj);
  NonLinSolver.setTolerance(1e-6);

  FomObserver<state_t> Obs("linadv_1d_solution.bin", 1);

  const auto dt = 0.001;
  const auto Nsteps = 2000;
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0., dt, Nsteps, Obs, NonLinSolver);

  return 0;
}
