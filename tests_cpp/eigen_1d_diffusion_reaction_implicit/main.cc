
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/diffusion_reaction.hpp"
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
  using state_t	= typename app_t::state_type;
  using rhs_t	= typename app_t::velocity_type;
  using jacob_t	= typename app_t::jacobian_type;

  auto state = appObj.initialCondition();
  auto stepperObj = pressio::ode::create_implicit_stepper(
    pressio::ode::StepScheme::BDF1, state, appObj);

  using lin_solver_t = pressio::linearsolvers::Solver<
    pressio::linearsolvers::iterative::Bicgstab, jacob_t>;
  lin_solver_t linSolverObj;
  auto NonLinSolver=
   pressio::nonlinearsolvers::create_newton_raphson(stepperObj, state, linSolverObj);
  NonLinSolver.setTolerance(1e-11);

  FomObserver<state_t> Obs("1d_diffreac_solution.bin", 10);

  const auto dt = 0.01;
  const auto Nsteps = 100;
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0., dt, Nsteps, Obs, NonLinSolver);

  return 0;
}
