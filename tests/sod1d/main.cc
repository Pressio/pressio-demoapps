
#include "pressio_ode_explicit.hpp"
#include "euler.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
{
  // pressio::log::initialize(pressio::logto::terminal);
  // pressio::log::setVerbosity({pressio::log::level::debug});

  using scalar_t = double;
  using mesh_t = pressiodemoapps::CellCenteredUniformMeshEigen<scalar_t>;
  using app_t       = pressiodemoapps::Sod1dEigen<scalar_t, mesh_t>;
  using app_state_t = typename app_t::state_type;
  using app_rhs_t   = typename app_t::velocity_type;

  mesh_t meshObj(".");

  app_t appObj(meshObj);
  const auto gamma = appObj.gamma();

  const auto stateSize = appObj.totalDofStencilMesh();
  using ode_state_t = pressio::containers::Vector<app_state_t>;
  ode_state_t state(appObj.initialCondition());

  auto stepperObj = pressio::ode::createRungeKutta4Stepper(state, appObj);
  FomObserver<ode_state_t> Obs("solution.bin", 1);

  const auto dt = 0.001;
  const auto Nsteps = 100;

  pressio::ode::advanceNSteps(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
