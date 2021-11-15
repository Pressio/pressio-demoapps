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

  const auto probId  = pda::Euler2d::SedovSymmetry;
  auto appObj      = pda::createExplicitProblemEigen(meshObj, probId, order);
  const auto stateSize = appObj.totalDofStencilMesh();

  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;

  state_t state = appObj.initialCondition();
  FomObserver<state_t> Obs("sedov2dsym_solution.bin", 10);
  auto stepperObj = pressio::ode::create_ssprk3_stepper(state, appObj);

  const auto dt = 0.0001;
  const auto Nsteps = 0.05/dt;
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
