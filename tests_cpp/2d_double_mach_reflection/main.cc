
#include "pressio_ode_explicit.hpp"
#include "euler2d.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
{

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
#if defined USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#else
  const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probId  = pda::Euler2d::DoubleMachReflection;
  auto appObj      = pda::createProblemEigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using app_state_t = typename app_t::state_type;
  using ode_state_t = pressio::containers::Vector<app_state_t>;

  const auto stateSize = appObj.totalDofStencilMesh();
  auto ic = appObj.initialCondition();
  ode_state_t state(ic);

  const auto dt = 0.0005;
  const auto Nsteps = 0.25/dt;
  auto time = 0.0;
  FomObserver<ode_state_t> Obs("doubleMach2d_solution.bin", 50);
  app_state_t state1(stateSize);
  app_state_t state2(stateSize);
  app_state_t veloTmp(stateSize);
  Obs(0, time, state);
  for (int step=1; step<=Nsteps; ++step)
  {
    if( step % 50 == 0){
      std::cout << "step = " << step << " out of " << Nsteps << std::endl;
    }

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
