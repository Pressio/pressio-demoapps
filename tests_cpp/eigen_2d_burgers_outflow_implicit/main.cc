
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/advection_diffusion2d.hpp"
#include "../observer.hpp"

int main()
{

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
#ifdef USE_WENO5
  const auto scheme   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  const auto scheme   = pda::InviscidFluxReconstruction::Weno3;
#else
  const auto scheme   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto schemeVisc = pda::ViscousFluxReconstruction::FirstOrder;

  const auto probId  = pda::AdvectionDiffusion2d::BurgersOutflow;
  auto appObj      = pda::create_problem_eigen(meshObj, probId, scheme, schemeVisc);

  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;
  using jacob_t	= typename app_t::jacobian_type;
  state_t state = appObj.initialCondition();

  auto stepperObj = pressio::ode::create_implicit_stepper(
    pressio::ode::StepScheme::CrankNicolson, appObj);

  using lin_solver_t = pressio::linearsolvers::Solver<
    pressio::linearsolvers::iterative::Bicgstab, jacob_t>;
  lin_solver_t linSolverObj;
  auto NonLinSolver=pressio::create_newton_solver(stepperObj, linSolverObj);
  NonLinSolver.setStopTolerance(1e-5);

  const auto dt = 0.01;
  const auto Nsteps = pressio::ode::StepCount(2./dt);
  FomObserver<state_t> Obs("burgers2d_solution.bin", 50);
  pressio::ode::advance_n_steps(stepperObj, state, 0., dt, Nsteps, Obs, NonLinSolver);

  return 0;
}
