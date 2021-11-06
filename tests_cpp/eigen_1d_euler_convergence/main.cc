
#include "pressio_ode_explicit.hpp"
#include "../observer.hpp"
#include "euler1d.hpp"

template<class scalar_type>
void analytical(const scalar_type x,
     scalar_type t,
     std::array<scalar_type,3> & prim)
{
  prim[0] = 1. + 0.2*std::sin(M_PI*(x-t) );
  prim[1] = 1.;
  prim[2] = 1.;
}

int main(int argc, char *argv[])
{

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");

#ifdef USE_WENO5
  const auto order   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#endif

  const auto probId  = pda::Euler1d::PeriodicSmooth;
  auto appObj	     = pda::createProblemEigen(meshObj, probId, order);
  // by default, velocoity and Jac are fused, but we
  // only want velocity here
  appObj.disableJacobian();

  using app_t = decltype(appObj);
  using app_state_t = typename app_t::state_type;
  using ode_state_t = pressio::containers::Vector<app_state_t>;

  using ode_state_t = pressio::containers::Vector<app_state_t>;
  ode_state_t state(appObj.initialCondition());

  auto stepperObj = pressio::ode::createRungeKutta4Stepper(state, appObj);
#ifdef USE_WENO5
  FomObserver<ode_state_t> Obs("eigen_1d_euler_convergence_weno5_solution.bin", 1);
#elif defined USE_WENO3
  FomObserver<ode_state_t> Obs("eigen_1d_euler_convergence_weno3_solution.bin", 1);
#endif

  const auto dt = 0.001;
  const auto Nsteps = 2./dt;
  pressio::ode::advanceNSteps(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
