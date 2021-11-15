
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
{

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
#ifdef USE_WENO5
  const auto order   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#else
  const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probId  = pda::Euler2d::SedovFull;
  auto appObj      = pda::createExplicitProblemEigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;

  const auto stateSize = appObj.totalDofStencilMesh();
  state_t state = appObj.initialCondition();

  FomObserver<state_t> Obs("sedov2d_solution.bin", 5);

  const auto dt = 0.0001;
  const auto Nsteps = 0.05/dt;
  auto stepperObj = pressio::ode::create_ssprk3_stepper(state, appObj);
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
