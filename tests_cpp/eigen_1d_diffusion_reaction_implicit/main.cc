
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/diffusion_reaction1d.hpp"
#include "../observer.hpp"

int main()
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
  auto appObj        = pda::create_diffusion_reaction_1d_problem_A_eigen(meshObj, 0.01, 0.005);
  using app_t = decltype(appObj);
  using state_t	= typename app_t::state_type;
  using jacob_t	= typename app_t::jacobian_type;

  auto state = appObj.initialCondition();
  auto stepperObj = pressio::ode::create_implicit_stepper(
    pressio::ode::StepScheme::CrankNicolson, appObj);

  using lin_solver_t = pressio::linearsolvers::Solver<
    pressio::linearsolvers::iterative::Bicgstab, jacob_t>;
  lin_solver_t linSolverObj;
  auto NonLinSolver=pressio::create_newton_solver(stepperObj, linSolverObj);
  NonLinSolver.setStopTolerance(1e-11);

  FomObserver<state_t> Obs("1d_diffreac_solution.bin", 10);

  const auto dt = 0.005;
  const auto Nsteps = pressio::ode::StepCount(200);
  pressio::ode::advance_n_steps(stepperObj, state, 0., dt, Nsteps, Obs, NonLinSolver);

  return 0;
}
