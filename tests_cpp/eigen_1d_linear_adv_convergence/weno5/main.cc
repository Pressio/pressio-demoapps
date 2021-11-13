
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/advection.hpp"
#include "../../observer.hpp"

// https://arxiv.org/pdf/1609.07625.pdf
// see page 12, example 1

template<class scalar_type>
scalar_type analytical(const scalar_type x, scalar_type t)
{
  const auto arg = M_PI*(x-t);
  return std::sin( arg );
}

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
  const auto order   = pda::InviscidFluxReconstruction::Weno5;
  const auto probid = pda::Advection1d::PeriodicLinear;
  auto appObj      = pda::createExplicitProblemEigen(meshObj, probid, order);

  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;

  const auto stateSize = appObj.totalDofStencilMesh();
  state_t state(stateSize);
  const auto & x = meshObj.viewX();
  for (int i=0; i<state.size(); ++i){
    state(i) = analytical(x(i), 0.0);
  }

  auto stepperObj = pressio::ode::create_rk4_stepper(state, appObj);
  FomObserver<state_t> Obs("eigen_1d_linear_adv_convergence_weno5_solution.bin", 10);

  const auto dt = 0.001;
  const auto Nsteps = 2./dt;
  pressio::ode::advance_n_steps_and_observe(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
