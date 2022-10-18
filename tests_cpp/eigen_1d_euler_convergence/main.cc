
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler1d.hpp"
#include "../observer.hpp"

template<class scalar_type>
void analytical(const scalar_type x,
     scalar_type t,
     std::array<scalar_type,3> & prim)
{
  prim[0] = 1. + 0.2*std::sin(M_PI*(x-t) );
  prim[1] = 1.;
  prim[2] = 1.;
}

int main()
{

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");

#ifdef USE_WENO5
  const auto order   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#endif

  const auto probId  = pda::Euler1d::PeriodicSmooth;
  auto appObj	     = pda::create_problem_eigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;
  state_t state(appObj.initialCondition());

  auto stepperObj = pressio::ode::create_rk4_stepper(appObj);
#ifdef USE_WENO5
  FomObserver<state_t> Obs("eigen_1d_euler_convergence_weno5_solution.bin", 1);
#elif defined USE_WENO3
  FomObserver<state_t> Obs("eigen_1d_euler_convergence_weno3_solution.bin", 1);
#endif

  const auto dt = 0.001;
  const auto Nsteps = pressio::ode::StepCount(2./dt);
  pressio::ode::advance_n_steps(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
