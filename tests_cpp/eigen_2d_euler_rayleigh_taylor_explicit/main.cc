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

  const double amplitude = 0.025;

  const auto probId  = pda::Euler2d::RayleighTaylor;
  auto appObj = pda::create_problem_eigen(meshObj, probId, order, amplitude);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;

  if (appObj.gamma() != 5./3.){
    throw std::runtime_error("RayleighTaylor: test: gamma is not correct");
  }

  state_t state(appObj.initialCondition());

  FomObserver<state_t> Obs("rt2d_solution.bin", 50);
#ifdef USE_WENO5
  const auto dt = 0.001;
#elif defined USE_WENO3
  const auto dt = 0.001;
#else
  const auto dt = 0.0005;
#endif

  std::cout << "dt = " << dt << "\n";
  const auto Nsteps = 1.95/dt;

  auto stepperObj = pressio::ode::create_ssprk3_stepper(state, appObj);
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0., dt, Nsteps, Obs);

  pressio::log::finalize();
  return 0;
}
