
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler3d.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
{

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
#ifdef USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#else
  const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probId  = pda::Euler3d::SedovSymmetry;
  auto appObj      = pda::create_problem_eigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;

  const auto stateSize = appObj.totalDofStencilMesh();
  state_t state(appObj.initialCondition());

  auto stepperObj = pressio::ode::create_ssprk3_stepper(state, appObj);
  FomObserver<state_t> Obs("3d_sedov_sym_solution.bin", 1);

  const auto dt = 0.0001;
  const auto Nsteps = 0.025/dt;
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
