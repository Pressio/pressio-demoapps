
#include "pressio_ode_explicit.hpp"
#include "euler2d.hpp"
#include "../observer.hpp"
#include <random>

std::mt19937 gen(2195884);
std::uniform_int_distribution<double> distrib(1.1, 2.0);

template<class scalar_type, class prim_t>
void initcond(int i, const scalar_type x, const scalar_type y, prim_t & prim)
{
  prim[0] = 1. + 0.1*x + 0.2*y;
  prim[1] = 0.001;
  prim[2] = 0.002;
  prim[3] = 0.00001;
}

template<class T>
void writeToFile(const T& obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (size_t i=0; i<obj.size(); i++){
    file << std::setprecision(14) << obj(i) << " \n";
  }
  file.close();
}

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
#ifdef USE_WENO5
  const auto order   = pda::ReconstructionType::fifthOrderWeno;
#else
  const auto order   = pda::ReconstructionType::firstOrder;
#endif

  const auto probId  = pda::Euler2d::testingonlyneumann;
  auto appObj      = pda::createProblemEigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using app_state_t = typename app_t::state_type;
  using ode_state_t = pressio::containers::Vector<app_state_t>;
  using scalar_t = typename app_t::scalar_type;

  const auto gamma = appObj.gamma();
  const auto stateSize = appObj.totalDofStencilMesh();
  app_state_t state(stateSize);
  std::array<scalar_t, 4> prim;
  const auto &x = meshObj.viewX();
  const auto &y = meshObj.viewY();
  for (int i=0; i<x.size(); ++i){
    initcond(i, x(i), y(i), prim);
    const auto ind = i*4;
    state(ind)   = prim[0];
    state(ind+1) = prim[0]*prim[1];
    state(ind+2) = prim[0]*prim[2];
    state(ind+3) = pressiodemoapps::ee::computeEnergyFromPrimitive(gamma, prim);
  }

  auto time = 0.0;
  auto velo = appObj.createVelocity();
  appObj.velocity(state, time, velo);
  writeToFile(velo,  "velo.txt");
  writeToFile(state, "state.txt");

  // for (int i=0; i<velo.size(); ++i){
  //   std::cout << i << " " << i % 4 << " " << std::setprecision(14) << velo(i) << "\n";
  // }

  return 0;
}
