
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/diffusion_reaction2d.hpp"
#include "../observer.hpp"

int main()
{

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
  const auto probId  = pda::DiffusionReaction2d::ProblemA;
  const auto recEn   = pda::ViscousFluxReconstruction::FirstOrder;
  auto appObj        = pda::create_problem_eigen(meshObj, probId, recEn);

  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;
  state_t state = appObj.initialCondition();

  const auto dt = 0.001;
  const auto Nsteps = pressio::ode::StepCount(100);
  FomObserver<state_t> Obs("diffusion_reaction_2d_solution.bin", 10);

  auto stepperObj = pressio::ode::create_ssprk3_stepper(appObj);
  pressio::ode::advance_n_steps(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
