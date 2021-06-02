
#include "pressio_ode_explicit.hpp"
#include "euler.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
{
  // pressio::log::initialize(pressio::logto::terminal);
  // pressio::log::setVerbosity({pressio::log::level::debug});

  using scalar_t = double;
  using mesh_t = pressiodemoapps::CellCenteredUniformMeshEigen<scalar_t>;
  using app_t       = pressiodemoapps::Riemann2dEigen<scalar_t, mesh_t>;
  using app_state_t = typename app_t::state_type;
  using app_rhs_t   = typename app_t::velocity_type;

  mesh_t meshObj(".");
  app_t appObj(meshObj, 2);
  const auto gamma = appObj.gamma();

  const auto stateSize = appObj.totalDofStencilMesh();
  using ode_state_t = pressio::containers::Vector<app_state_t>;
  auto ic = appObj.initialCondition();
  ode_state_t state(ic);

  auto stepperObj = pressio::ode::createRungeKutta4Stepper(state, appObj);
  FomObserver<ode_state_t> Obs("solution.bin", 100);

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
