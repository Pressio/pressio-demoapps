
#include "pressio_ode_explicit.hpp"
#include "euler2d.hpp"
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

  const auto probId  = pda::Euler2d::Riemann;
  auto appObj      = pda::createProblemEigen(meshObj, probId, order, 2);
  appObj.disableJacobian();
  using app_t = decltype(appObj);
  using app_state_t = typename app_t::state_type;
  using ode_state_t = pressio::containers::Vector<app_state_t>;

  const auto stateSize = appObj.totalDofStencilMesh();
  ode_state_t state(appObj.initialCondition());

  FomObserver<ode_state_t> Obs("riemann2d_solution.bin", 100);
  const auto dt = 0.001;
  const auto Nsteps = 0.6/dt;

  auto time = 0.0;
  app_state_t state1(stateSize);
  app_state_t state2(stateSize);
  Obs(0, time, state);
  for (int step=1; step<=Nsteps; ++step)
  {
    app_state_t veloTmp(stateSize);
    appObj.velocity(*state.data(), time, veloTmp);
    state1 = *state.data() + dt*veloTmp;

    appObj.velocity(state1, time, veloTmp);
    state2 = (3./4.)*(*state.data()) + 0.25*state1 + 0.25*dt*veloTmp;

    appObj.velocity(state2, time, veloTmp);
    *state.data() = (1./3.)*((*state.data()) + 2.*state2 + 2.*dt*veloTmp);

    time = static_cast<double>(step) * dt;
    Obs(step, time, state);
  }

  return 0;
}
