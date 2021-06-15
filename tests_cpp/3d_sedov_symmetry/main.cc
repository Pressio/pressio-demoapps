
#include "pressio_ode_explicit.hpp"
#include "euler3d.hpp"
#include "../observer.hpp"

// template<class T>
// void writeToFile(const T& obj, const std::string & fileName)
// {
//   std::ofstream file; file.open(fileName);
//   for (size_t i=0; i<obj.extent(0); i++){
//     file << std::setprecision(14) << obj(i) << " \n";
//   }
//   file.close();
// }

int main(int argc, char *argv[])
{

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
  const auto order   = pda::ReconstructionType::firstOrder;

  const auto probId  = pda::Euler3d::SedovSymmetry;
  auto appObj      = pda::create3dProblemEigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using app_state_t = typename app_t::state_type;
  using ode_state_t = pressio::containers::Vector<app_state_t>;
  using scalar_t = typename app_t::scalar_type;

  const auto stateSize = appObj.totalDofStencilMesh();
  ode_state_t state(appObj.initialCondition());
  // writeToFile(state,  "ic.txt");

  auto stepperObj = pressio::ode::createRungeKutta4Stepper(state, appObj);
  FomObserver<ode_state_t> Obs("3d_sedov_sym_solution.bin", 1);

  const auto dt = 0.0001;
  const auto Nsteps = 0.025/dt;

  auto time = 0.0;
  app_state_t state1(stateSize);
  app_state_t state2(stateSize);
  app_state_t veloTmp(stateSize);
  Obs(0, time, state);
  for (int step=1; step<=Nsteps; ++step)
  {
    if( step % 50 == 0){
      std::cout << "step = " << step << std::endl;
    }

    appObj.velocity(*state.data(), time, veloTmp);
    state1 = *state.data() + dt*veloTmp;

    appObj.velocity(state1, time, veloTmp);
    state2 = (3./4.)*(*state.data()) + 0.25*state1 + 0.25*dt*veloTmp;

    appObj.velocity(state2, time, veloTmp);
    *state.data() = (1./3.)*((*state.data()) + 2.*state2 + 2.*dt*veloTmp);

    time = static_cast<double>(step) * dt;
    Obs(step, time, state);
  }

  return 0;
}
