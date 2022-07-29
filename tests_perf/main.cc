
#include "pressiodemoapps/swe2d.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "pressiodemoapps/euler3d.hpp"
#include <chrono>
#include "CLI11.hpp"

namespace pda = pressiodemoapps;

pda::InviscidFluxReconstruction stringToInviscidRecScheme(const std::string & string)
{
  if (string == "FirstOrder"){
    return pda::InviscidFluxReconstruction::FirstOrder;
  }
  else if (string == "Weno3"){
    return pda::InviscidFluxReconstruction::Weno3;
  }
  else if (string == "Weno5"){
    return pda::InviscidFluxReconstruction::Weno5;
  }

  return {};
}

int main(int argc, char *argv[])
{
  CLI::App app;

  std::string meshDir = "void";
  int loopCount = 10;
  std::string scheme = "void";
  std::string problem = "void";

  app.add_option("-m, --meshDir", meshDir);
  app.add_option("-n, --trialsCount", loopCount);
  app.add_option("-s, --scheme", scheme);
  //app.add_option("-p, --problem", problem);
  CLI11_PARSE(app, argc, argv);

  auto probE   = pda::Euler2d::PeriodicSmooth;
  auto schemeE = stringToInviscidRecScheme(scheme);

  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(meshDir);
  auto appObj = pda::create_problem_eigen(meshObj, probE, schemeE);
  using app_t = decltype(appObj);

  auto state = appObj.initialCondition();
  auto V     = appObj.createRightHandSide();
  auto J     = appObj.createJacobian();

  // warmup
  // appObj.rightHandSideAndJacobian(state, 0., V, J);
  appObj.rightHandSide(state, 0., V);

  auto t1 = std::chrono::high_resolution_clock::now();
  for (int i=0; i<loopCount; ++i){
    appObj.rightHandSide(state, 0., V);
    // appObj.rightHandSideAndJacobian(state, 0., V, J);
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration< double > fs = t2 - t1;
  std::cout << "elapsed " << fs.count() << std::endl;

  return 0;
}
