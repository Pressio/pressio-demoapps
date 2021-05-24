
#include "pressio_ode_explicit.hpp"
#include "advection.hpp"
#include "../observer.hpp"

// https://arxiv.org/pdf/1609.07625.pdf
// see page 12, example 1
// 

template<class scalar_type>
scalar_type analytical(const scalar_type x, scalar_type t)
{
  // const auto arg = M_PI*(x-t);
  // return std::sin( arg - (1./M_PI)*std::sin(arg) );
  const auto arg = M_PI*(x-t);
  return std::sin( arg );
}

int main(int argc, char *argv[])
{
  // pressio::log::initialize(pressio::logto::terminal);
  // pressio::log::setVerbosity({pressio::log::level::debug});

  using scalar_t = double;
  using mesh_t = pressiodemoapps::CellCenteredUniformMeshEigen<scalar_t>;
  using app_t       = pressiodemoapps::PeriodicLinearAdvection1dEigen<scalar_t, mesh_t>;
  using app_state_t = typename app_t::state_type;
  using app_rhs_t   = typename app_t::velocity_type;

  mesh_t meshObj(".");
  app_t appObj(meshObj);

  const auto stateSize = appObj.totalDofStencilMesh();
  using ode_state_t = pressio::containers::Vector<app_state_t>;
  ode_state_t state(stateSize);
  const auto & x = meshObj.viewX();
  for (int i=0; i<state.extent(0); ++i){
    state(i) = analytical(x(i), 0.0);
  }

  auto stepperObj = pressio::ode::createRungeKutta4Stepper(state, appObj);
  FomObserver<ode_state_t> Obs("solution.bin", 10);

  const auto dt = 0.001;
  const auto Nsteps = 2./dt;
  pressio::ode::advanceNSteps(stepperObj, state, 0., dt, Nsteps, Obs);

  // auto rhsStr = testRhs(appObj);
  // if (rhsStr){
  //   std::puts("PASS");
  // }
  // else{
  //   std::puts("FAIL");
  // }

  return 0;
}




// template<class T>
// bool testRhs(T & appObj)
// {
//   using scalar_type = typename T::scalar_type;

//   typename T::velocity_type fomVelo(appObj.createVelocity());
//   const auto & x = appObj.viewX();
//   const auto & y = appObj.viewY();
//   std::cout << x.size() << " " << y.size() << std::endl;
//   typename T::state_type fomState(fomVelo.size());
//   std::cout << fomState.size() << std::endl;

//   std::array<scalar_type, 4> prim;
//   for (int i=0; i<x.size(); ++i){
//     analytical(x(i), y(i), 0.0, prim);
//     const auto ind = i*4;
//     fomState(ind)   = prim[0];
//     fomState(ind+1) = prim[0]*prim[1];
//     fomState(ind+2) = prim[0]*prim[2];
//     fomState(ind+3) = computeEnergy(prim);
//   }

//   std::ofstream file; file.open("state.txt");
//   for (int i=0; i<fomState.size(); i++){
//     file << std::setprecision(14) << fomState(i) << " \n";
//   }
//   file.close();


//   appObj.velocity(fomState, 0., fomVelo);
//   std::cout << fomVelo.size() << std::endl;
//   std::cout << std::setprecision(14) << fomVelo << std::endl;

//   {
//   std::ofstream file; file.open("velo.txt");
//   for (int i=0; i<fomVelo.size(); i++){
//     file << std::setprecision(14) << fomVelo(i) << " \n";
//   }
//   file.close();
//   }

//   // typename T::state_type goldVelo(fomVelo.size());
//   // goldVelo <<
//   //   440.19802377938, 460.99843957979, 739.39885538021,
//   //   356.61259942078, 346.77688513507, 696.24117084935,
//   //   369.02281312333, 349.10531006335, 742.18780700338,
//   //   344.18175264043, 355.57260526606, 749.6634578917,
//   //   48.461276134967, 110.76718522588, -98.526905683215,
//   //   -44.340216773232, -89.634490920138, -232.82876506704,
//   //   -44.468414313835, -98.89192709138, -269.51543986893,
//   //   -98.573684863215, -47.93991491799, -359.80614497276,
//   //   66.592185379923, 162.47916453749, -230.43385630494,
//   //   -46.319478775087, -127.15399490412, -379.08851103315,
//   //   -46.880336994272, -136.45022889228, -415.42012079029,
//   //   -114.93204009409, -27.252662304976, -519.27328451586,
//   //   -183.7999939098, -54.746515648928, -31.693037388059,
//   //   -310.94451610532, -427.59411653631, -116.54371696731,
//   //   -321.93128757954, -446.86716215377, -90.403036728008,
//   //   -412.87868107042, -288.36857643775, -144.75847180507;

//   // if (!goldVelo.isApprox(fomVelo)){
//   //return false;
//   // }
//   // else{
//    return true;
//   // }
// }
