
#include "pressio_ode_explicit.hpp"
#include "diffusion_reaction.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
  const auto sch     = ::pressiodemoapps::ViscousFluxReconstruction::FirstOrder;
  auto appObj        = pda::createProblemEigen(meshObj,
					       pda::DiffusionReaction1d::ProblemA,
					       sch, 0.01, 0.005);
  using app_t = decltype(appObj);
  using app_state_t = typename app_t::state_type;
  using ode_state_t = pressio::containers::Vector<app_state_t>;
  ode_state_t state = appObj.initialCondition();

  auto stepperObj = pressio::ode::createRungeKutta4Stepper(state, appObj);
  FomObserver<ode_state_t> Obs("1d_diffreac_solution.bin", 10);

  const auto dt = 0.001;
  const auto Nsteps = 1000;
  pressio::ode::advanceNSteps(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
