
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/advection_diffusion2d.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
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

  const auto probId  = pda::AdvectionDiffusion2d::BurgersPeriodic;
  auto appObj      = pda::create_problem_eigen(meshObj, probId, scheme, schemeVisc);
  const auto stateSize = appObj.totalDofStencilMesh();

  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;
  state_t state = appObj.initialCondition();

  const auto dt = 0.005;
  const auto Nsteps = 10./dt;
  FomObserver<state_t> Obs("burgers2d_solution.bin", 50);

  auto stepperObj = pressio::ode::create_rk4_stepper(state, appObj);
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
