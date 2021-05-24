
#include "pressio_ode_explicit.hpp"
#include "euler.hpp"
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

int main(int argc, char *argv[])
{
  // pressio::log::initialize(pressio::logto::terminal);
  // pressio::log::setVerbosity({pressio::log::level::debug});

  using scalar_t = double;
  using mesh_t = pressiodemoapps::CellCenteredUniformMeshEigen<scalar_t>;
  using app_t       = pressiodemoapps::PeriodicEuler1dEigen<scalar_t, mesh_t>;
  using app_state_t = typename app_t::state_type;
  using app_rhs_t   = typename app_t::velocity_type;

  mesh_t meshObj(".");

  app_t appObj(meshObj);
  const auto gamma = appObj.gamma();

  const auto stateSize = appObj.totalDofStencilMesh();
  using ode_state_t = pressio::containers::Vector<app_state_t>;
  ode_state_t state(stateSize);
  std::array<scalar_t, 3> prim;
  const auto x = meshObj.viewX();
  for (int i=0; i<x.size(); ++i){
    analytical(x(i), 0.0, prim);
    const auto ind = i*3;
    state(ind)   = prim[0];
    state(ind+1) = prim[0]*prim[1];
    state(ind+2) = pressiodemoapps::ee::computeEnergyFromPrimitive(gamma, prim);
  }

  auto stepperObj = pressio::ode::createRungeKutta4Stepper(state, appObj);
  FomObserver<ode_state_t> Obs("solution.bin", 1);

  const auto dt = 0.001;
  const auto Nsteps = 2./dt;
  pressio::ode::advanceNSteps(stepperObj, state, 0., dt, Nsteps, Obs);

  return 0;
}
