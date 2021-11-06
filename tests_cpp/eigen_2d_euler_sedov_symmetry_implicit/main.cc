
#include "pressio_ode_implicit.hpp"
#include "euler2d.hpp"
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

  const auto probId  = pda::Euler2d::SedovSymmetry;
  auto appObj = pda::createProblemEigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  using app_jacob_t	= typename app_t::jacobian_type;

  using ode_state_t = pressio::containers::Vector<app_state_t>;
  using ode_res_t   = pressio::containers::Vector<app_rhs_t>;
  using ode_jac_t   = pressio::containers::SparseMatrix<app_jacob_t>;

  ode_state_t state = appObj.initialCondition();

  using ode_tag = pressio::ode::implicitmethods::CrankNicolson;
  using stepper_t = pressio::ode::ImplicitStepper<
    ode_tag, ode_state_t, ode_res_t, ode_jac_t, app_t>;
  stepper_t stepperObj(state, appObj);

  using lin_solver_t = pressio::solvers::linear::Solver<
    pressio::solvers::linear::iterative::Bicgstab, ode_jac_t>;
  lin_solver_t linSolverObj;
  auto NonLinSolver=
    pressio::solvers::nonlinear::createNewtonRaphson(stepperObj, state, linSolverObj);
  NonLinSolver.setTolerance(1e-5);

  FomObserver<ode_state_t> Obs("sedov2dsym_solution.bin", 5);

  const auto dt = 0.0001;
  const auto Nsteps = 100;
  pressio::ode::advanceNSteps(stepperObj, state, 0., dt, Nsteps, Obs, NonLinSolver);

  return 0;
}
