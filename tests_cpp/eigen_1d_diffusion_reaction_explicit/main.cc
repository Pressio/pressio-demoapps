
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/diffusion_reaction.hpp"
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
  using state_t = typename app_t::state_type;
  auto state = appObj.initialCondition();

  auto stepperObj = pressio::ode::create_rk4_stepper(state, appObj);
  FomObserver<state_t> Obs("1d_diffreac_solution.bin", 10);

  const auto dt = 0.001;
  const auto Nsteps = 1000;
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
