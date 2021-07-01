
#include "pressio_ode_explicit.hpp"
#include "euler2d.hpp"
#include "observer.hpp"
#include "CLI11.hpp"

int main(int argc, char *argv[])
{
   CLI::App app{"2D Double Mach reflection- FOM"};

  using scalar_t = double;
  scalar_t dt = 0.0005;
  int sampleEvery = 50;
  scalar_t finalTime = 0.25;
  app.add_option("--sampleFreq", sampleEvery,
		 "Frequency of sampling state");

  app.add_option("--dt", dt,
		 "Time step size: default = 0.0001");

  app.add_option("-T,--finalTime", finalTime,
		 "Simulation time: default = 0.25");
  CLI11_PARSE(app, argc, argv);

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
  const auto probId  = pda::Euler2d::DoubleMachReflection;
  auto appObj      = pda::createProblemEigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using app_state_t = typename app_t::state_type;
  using ode_state_t = pressio::containers::Vector<app_state_t>;
  using scalar_t = typename app_t::scalar_type;

  const auto stateSize = appObj.totalDofStencilMesh();
  auto ic = appObj.initialCondition();
  ode_state_t state(ic);

  FomObserver<ode_state_t> Obs("fom_dmr_solution.bin", sampleEvery);
  const auto Nsteps = finalTime/dt;
  auto time = 0.0;
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
