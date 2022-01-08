
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/diffusion_reaction.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  namespace pda    = pressiodemoapps;
  namespace plsol  = pressio::linearsolvers;
  namespace pnlsol = pressio::nonlinearsolvers;

  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
  const auto scheme  = pda::ViscousFluxReconstruction::FirstOrder;
  const auto probId  = pda::DiffusionReaction3d::GrayScott;
  auto appObj        = pda::create_problem_eigen(meshObj, probId, scheme);
  const auto stateSize = appObj.totalDofStencilMesh();

  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;
  using jacob_t = typename app_t::jacobian_type;
  state_t state = appObj.initialCondition();

  auto stepperObj = pressio::ode::create_bdf1_stepper(state, appObj);

  using lin_solver_t = plsol::Solver<plsol::iterative::Bicgstab, jacob_t>;
  lin_solver_t linSolverObj;
  auto NonLinSolver= pnlsol::create_newton_raphson(stepperObj, state, linSolverObj);
  NonLinSolver.setTolerance(1e-5);

  const auto dt = 1.;
  const std::size_t Nsteps = 100./dt;
  FomObserver<state_t> Obs("gs_3d_solution.bin", 10);
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0., dt, Nsteps, Obs, NonLinSolver);

  pressio::log::finalize();
  return 0;
}
