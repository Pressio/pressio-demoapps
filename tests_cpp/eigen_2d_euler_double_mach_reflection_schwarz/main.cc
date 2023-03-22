
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "../observer.hpp"

int main()
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");

  constexpr auto order   = pda::InviscidFluxReconstruction::FirstOrder;
  const auto probId  = pda::Euler2d::DoubleMachReflection;

  auto appObj = pda::create_problem_eigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;
  using jacob_t = typename app_t::jacobian_type;

  state_t state = appObj.initialCondition();

  int numDomains = 4;
  int nx = 50;
  int ny = 50;
  state_t v((nx + ny ) * 2 * 4);
  std::vector<state_t> vecStateBcs( numDomains, v);

  appObj.setStateBc(&vecStateBcs[0]);
  // for (int domIdx = 0; domIdx < numDomains; ++domIdx){
  //   vecApp[domIdx].setStateBc(vecStateBcs[domIdx]);
  // }

  auto stepperObj = pressio::ode::create_implicit_stepper(
    pressio::ode::StepScheme::CrankNicolson, appObj);
  using lin_solver_t = pressio::linearsolvers::Solver<
    pressio::linearsolvers::iterative::Bicgstab, jacob_t>;
  lin_solver_t linSolverObj;
  auto NonLinSolver = pressio::nonlinearsolvers::create_newton_raphson(stepperObj, linSolverObj);
  NonLinSolver.setTolerance(1e-5);

  FomObserver<state_t> Obs("doubleMach2d_solution.bin", 10);

  // for (int n = 0; n < numSteps; ++n){
  //     double t = n * dt;
  //     while (conv > tol){
  // 	  for (int domIdx = 0; domIdx < numDomains; ++domIdx){
  // 	    advance_n_steps(stepperVec[domIdx], stateVec[domIdx], t, dt,
  // 			    1, Obs,	solverVec[domIdx]);
  // 	    // update stateBc using stateVec
  // 	  }
  // 	  conv = calcConvergence(stateVec);
  // 	}
  //   }

  const auto dt = 0.0025;
  const auto Nsteps = pressio::ode::StepCount(100);
  pressio::ode::advance_n_steps(stepperObj, state, 0., dt, Nsteps, Obs, NonLinSolver);

  return 0;
}
