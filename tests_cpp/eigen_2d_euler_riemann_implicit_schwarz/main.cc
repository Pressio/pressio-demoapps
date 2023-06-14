
#include "pressiodemoapps/euler2d.hpp"
#include "pressiodemoapps/impl/schwarz.hpp"
#include "../observer.hpp"

using namespace std;

int main()
{
  namespace pda  = pressiodemoapps;
  namespace plog = pressio::log;
  namespace pode = pressio::ode;

  plog::initialize(pressio::logto::terminal);
  plog::setVerbosity({pressio::log::level::debug});

  // +++++ USER INPUTS +++++
  string meshRoot = "/home/crwentl/research/runs/pressio/riemann/meshes/mesh_2x2_refx";
  string obsRoot = "riemann2d_solution";
  const int obsFreq = 1;

  // problem definition
  const auto probId = pda::Euler2d::Riemann;
  const auto order  = pda::InviscidFluxReconstruction::FirstOrder;
  const auto scheme = pode::StepScheme::CrankNicolson;

  // time stepping
  const int numSteps = 1000;
  vector<double> dt(1, 0.001);
  const int convergeStepMax = 10;
  const double abs_err_tol = 1e-11;
  const double rel_err_tol = 1e-11;

  // +++++ END USER INPUTS +++++

  // decomposition
  auto decomp = pda::impl::SchwarzDecomp<
    decltype(probId),
    decltype(order),
    decltype(scheme)
  >(probId, order, scheme, meshRoot, dt, 2);

  // observer
  using state_t = decltype(decomp)::state_t;
  using obs_t = FomObserver<state_t>;
  vector<obs_t> obsVec(decomp.ndomains);
  for (int domIdx = 0; domIdx < decomp.ndomains; ++domIdx) {
    obsVec[domIdx] = obs_t(obsRoot + "_" + to_string(domIdx) + ".bin", 1);
    obsVec[domIdx](::pressio::ode::StepCount(0), 0.0, decomp.stateVec[domIdx]);
  }

  // solve
  double time = 0.0;
  for (int outerStep = 1; outerStep <= numSteps; ++outerStep)
  {
    cout << "Step " << outerStep << endl;

    // compute contoller step until convergence
    decomp.calc_controller_step(
      outerStep,
      time,
      rel_err_tol,
      abs_err_tol,
      convergeStepMax
    );

    time += decomp.dtMax;

    // output observer
    const auto stepWrap = pode::StepCount(outerStep);
    for (int domIdx = 0; domIdx < decomp.ndomains; ++domIdx) {
      obsVec[domIdx](stepWrap, time, decomp.stateVec[domIdx]);
    }

  }

  return 0;
}
