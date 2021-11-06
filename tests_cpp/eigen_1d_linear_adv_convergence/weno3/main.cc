
#include "pressio_ode_explicit.hpp"
#include "advection.hpp"
#include "../../observer.hpp"

// https://link.springer.com/content/pdf/10.1007/s12591-019-00508-5.pdf

template<class scalar_type>
scalar_type analytical(const scalar_type x, scalar_type t)
{
  const auto pixmt = M_PI*(x-t);
  const auto sinval = std::sin(pixmt);
  return sinval*sinval*sinval*sinval;
}

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
  const auto probid = pda::Advection1d::PeriodicLinear;
  auto appObj      = pda::createProblemEigen(meshObj, probid, order);

  using app_t = decltype(appObj);
  using app_state_t = typename app_t::state_type;
  using ode_state_t = pressio::containers::Vector<app_state_t>;

  const auto stateSize = appObj.totalDofStencilMesh();
  ode_state_t state(stateSize);
  const auto & x = meshObj.viewX();
  for (int i=0; i<state.extent(0); ++i){
    state(i) = analytical(x(i), 0.0);
  }

  auto stepperObj = pressio::ode::createRungeKutta4Stepper(state, appObj);
  FomObserver<ode_state_t> Obs("eigen_1d_linear_adv_convergence_weno3_solution.bin", 10);

  const auto dt = 0.0001;
  const auto Nsteps = 2./dt;
  pressio::ode::advanceNSteps(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
