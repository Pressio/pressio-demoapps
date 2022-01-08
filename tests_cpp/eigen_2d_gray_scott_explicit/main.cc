
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/diffusion_reaction.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::info});

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
  const auto scheme  = pda::ViscousFluxReconstruction::FirstOrder;
  const auto probId  = pda::DiffusionReaction2d::GrayScott;
  auto appObj        = pda::create_problem_eigen(meshObj, probId, scheme);
  const auto stateSize = appObj.totalDofStencilMesh();

  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;
  state_t state = appObj.initialCondition();

  const auto dt = 0.8;
  const std::size_t Nsteps = 4000./dt;
  FomObserver<state_t> Obs("gs_2d_solution.bin", 20);

  auto stepperObj = pressio::ode::create_rk4_stepper(state, appObj);
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0., dt, Nsteps, Obs);

  pressio::log::finalize();
  return 0;
}
