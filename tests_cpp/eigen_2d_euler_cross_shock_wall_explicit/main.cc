
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
{

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
#ifdef USE_WENO5
  const auto recEnum   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  const auto recEnum   = pda::InviscidFluxReconstruction::Weno3;
#else
  const auto recEnum   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  auto appObj = pda::create_cross_shock_wall_problem_eigen(meshObj, recEnum, 0.05, 10.);
  const auto stateSize = appObj.totalDofStencilMesh();

  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;
  state_t state = appObj.initialCondition();

  auto time = 0.0;
  FomObserver<state_t> Obs("eulerCrossShock2dwall_solution.bin", 100);

  auto stepperObj = pressio::ode::create_rk4_stepper(state, appObj);

  const auto dt1 = 0.0005;
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0., dt1, 100./dt1, Obs);

  // const auto dt2 = 0.001;
  // pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0.05, dt2, 5.2/dt2, Obs);

  // const auto dt3 = 0.001;
  // pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0.25, dt3, 10./dt3, Obs);

  return 0;
}
