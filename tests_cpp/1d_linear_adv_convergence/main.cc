
#include "pressio_ode_explicit.hpp"
#include "advection.hpp"
#include "../observer.hpp"

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
  const auto order   = pda::reconstructionEnum::fifthOrderWeno;
  auto appObj	     = pda::createPeriodicLinearAdvection1dEigen(meshObj, order);

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
  FomObserver<ode_state_t> Obs("1d_linear_adv_convergence_solution.bin", 10);

  const auto dt = 0.001;
  const auto Nsteps = 2./dt;
  pressio::ode::advanceNSteps(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
