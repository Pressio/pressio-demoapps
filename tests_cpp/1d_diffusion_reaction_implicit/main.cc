
#include "pressio_ode_implicit.hpp"
#include "diffusion_reaction.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
  const auto sch     = ::pressiodemoapps::ViscousFluxReconstruction::FirstOrder;
  auto appObj        = pda::createProblemEigen(meshObj,
					       pda::DiffusionReaction1d::ProblemA,
					       sch, 0.01, 0.005);
  using app_t = decltype(appObj);

  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  using app_jacob_t	= typename app_t::jacobian_type;

  using ode_state_t = pressio::containers::Vector<app_state_t>;
  using ode_res_t   = pressio::containers::Vector<app_rhs_t>;
  using ode_jac_t   = pressio::containers::SparseMatrix<app_jacob_t>;

  ode_state_t state = appObj.initialCondition();

  using ode_tag = pressio::ode::implicitmethods::Euler;
  using stepper_t = pressio::ode::ImplicitStepper<
    ode_tag, ode_state_t, ode_res_t, ode_jac_t, app_t>;
  stepper_t stepperObj(state, appObj);

  using lin_solver_t = pressio::solvers::linear::Solver<
    pressio::solvers::linear::iterative::Bicgstab, ode_jac_t>;
  lin_solver_t linSolverObj;
  auto NonLinSolver=
    pressio::solvers::nonlinear::createNewtonRaphson(stepperObj, state, linSolverObj);
  NonLinSolver.setTolerance(1e-11);

  FomObserver<ode_state_t> Obs("1d_diffreac_solution.bin", 10);

  const auto dt = 0.01;
  const auto Nsteps = 100;
  pressio::ode::advanceNSteps(stepperObj, state, 0., dt, Nsteps, Obs, NonLinSolver);

  return 0;
}
