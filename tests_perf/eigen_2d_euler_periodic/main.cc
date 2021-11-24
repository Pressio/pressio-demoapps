
#include "pressiodemoapps/euler2d.hpp"
#include <chrono>

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(std::string(argv[1]));
#ifdef USE_WENO5
  const auto order   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#else
  const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  // here we want to test Euler2d with periodic BC, so we can use SmoothProblem
  // since that has periodic BC. THe problem itself does not really
  // matter as long as it is one such that it meets the target periodic BC.
  const auto probId = pda::Euler2d::PeriodicSmooth;
  auto appObj       = pda::create_problem_eigen(meshObj, probId, order);
  using app_t       = decltype(appObj);
  using state_t	    = typename app_t::state_type;
  state_t state(appObj.initialCondition());

  auto V = appObj.createVelocity();
  auto J = appObj.createJacobian();
  // warmup
  appObj.velocityAndJacobian(state, 0., V, J);
  //appObj.velocity(state, 0., V);

  const int N = std::stoi(argv[2]);
  auto t1 = std::chrono::high_resolution_clock::now();
  for (int i=0; i<N; ++i){
    appObj.velocityAndJacobian(state, 0., V, J);
    //appObj.velocity(state, 0., V);
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration< double > fs = t2 - t1;
  std::cout << "elapsed  = " << fs.count() << std::endl;
  std::cout << "per_iter = " << fs.count()/(double)N << std::endl;

  return 0;
}
