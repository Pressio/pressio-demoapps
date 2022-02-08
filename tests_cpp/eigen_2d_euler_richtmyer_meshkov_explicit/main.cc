#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "../observer.hpp"

int main(int argc, char *argv[])
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
#ifdef USE_WENO5
  const auto order   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#else
  const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const double gamma = 1.276;
  const double amplitude = 0.1;

  const auto probId  = pda::Euler2d::RichtmyerMeshkov;
  auto appObj      = pda::create_problem_eigen(meshObj, probId, order,
					       gamma, amplitude);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;

  state_t state(appObj.initialCondition());

  FomObserver<state_t> Obs("rmi2d_solution.bin", 5);
  const double cfl = 0.45;
  const auto dx = meshObj.dx();
  const auto dt = 0.002; //dx*cfl;
  std::cout << "dt = " << dt << "\n";
  const auto Nsteps = 8.0/dt;

  auto stepperObj = pressio::ode::create_ssprk3_stepper(state, appObj);
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0., dt, Nsteps, Obs);

  pressio::log::finalize();
  return 0;
}
