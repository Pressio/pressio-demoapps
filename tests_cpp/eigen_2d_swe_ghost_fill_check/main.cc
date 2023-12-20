#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/swe2d.hpp"

auto printLam = [](std::string s,  int it, auto const & M){
  std::cout << '\n';
  std::cout << s << ": ";
  for (int j=0; j< M.cols(); ++j){
    const auto v = M(it,j);

    if (v == std::numeric_limits<double>::min()){
      std::cout << "- ";
    }
    else{
      std::cout << v << " ";
    }
  }
 };

bool test_sample_s5()
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("./sample_mesh_s5");
  const auto order   = pda::InviscidFluxReconstruction::Weno3;

  const auto probId  = pda::Swe2d::SlipWall;
  auto appObj      = pda::create_problem_eigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;

  state_t state(appObj.initialCondition());
  double count = 0.;
  for (int i=0; i<state.size(); i+=3){
    state(i) = count;
    state(i+1) = count;
    state(i+2) = count++;
  }

  auto f = appObj.createRhs();
  appObj.rhs(state, 0., f);

  const auto & GL = appObj.viewGhostLeft();
  const auto & GF = appObj.viewGhostFront();
  const auto & GR = appObj.viewGhostRight();
  const auto & GB = appObj.viewGhostBack();

  auto checkLam = [](auto const & from, auto const & gold) ->bool{
    for (int i=0; i< gold.size(); ++i){
      if (gold[i] != from[i]){ return false; }
    }
    return true;
  };

  constexpr auto v = std::numeric_limits<double>::min();
  std::map<int, std::array<double, 6>> goldL;
  goldL[0]  = {0, -0, 0, 1, -1, 1};
  goldL[21] = {21, -21, 21, 22,-22,22};
  goldL[22] = {v,v,v, 21,-21,21};
  goldL[42] = {42,-42,42,43,-43,43};

  std::map<int, std::array<double, 6>> goldF;
  goldF[37] = {v,v,v,45,45,-45};
  goldF[42] = {42,42,-42,34,34,-34};
  goldF[45] = {45,45,-45,37,37,-37};
  goldF[49] = {49,49,-49,41,41,-41};

  std::map<int, std::array<double, 6>> goldR;
  goldR[7]  = {7,-7,7,6,-6,6};
  goldR[27] = {v,v,v,28,-28,28};
  goldR[28] = {28,-28,28,27,-27,27};
  goldR[49] = {49,-49,49,48,-48,48};

  std::map<int, std::array<double, 6>> goldB;
  goldB[0]  = {0,0,-0,8,8,-8};
  goldB[3]  = {3,3,-3,11,11,-11};
  goldB[7]  = {7,7,-7,15,15,-15};
  goldB[11] = {v,v,v,3,3,-3};

  const auto & graph = meshObj.graph();
  const auto & rowsBd = meshObj.graphRowsOfCellsNearBd();
  for (decltype(rowsBd.size()) it=0; it<rowsBd.size(); ++it){
    int smPt = rowsBd[it];
    const auto cGID = graph(smPt, 0);
    std::cout << cGID << ": ";

    printLam("left   ", it, GL);
    if (cGID == 0 || cGID==21 || cGID==22 || cGID==42){
      if (!checkLam( GL.row(it), goldL[cGID] )){
    	return false;
      }
    }

    printLam("back   ", it, GB);
    if (cGID == 0 || cGID==3 || cGID==7 || cGID==11){
      if (!checkLam( GB.row(it), goldB[cGID] )){
    	return false;
      }
    }

    printLam("front  ", it, GF);
    if (cGID == 37 || cGID==42 || cGID==45 || cGID==49){
      if (!checkLam( GF.row(it), goldF[cGID] )){
    	return false;
      }
    }

    printLam("right  ", it, GR);
    if (cGID == 7 || cGID==27 || cGID==28 || cGID==49){
      if (!checkLam( GR.row(it), goldR[cGID] )){
    	return false;
      }
    }

    std::cout << "\n";
    // std::cout << "\n";
  }

  return true;
}

// bool test_sample_s7()
// {
//   namespace pda = pressiodemoapps;
//   const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("./sample_mesh_s7");
//   const auto order   = pda::InviscidFluxReconstruction::Weno5;

//   const auto probId  = pda::Euler2d::Riemann;
//   auto appObj      = pda::create_problem_eigen(meshObj, probId, order, 1);
//   using app_t = decltype(appObj);
//   using state_t = typename app_t::state_type;

//   state_t state(appObj.initialCondition());
//   double count = 0.;
//   for (int i=0; i<state.size(); i+=4){
//     state(i) = count;
//     state(i+1) = count;
//     state(i+2) = count;
//     state(i+3) = count++;
//   }

//   auto f = appObj.createRhs();
//   appObj.rhs(state, 0., f);

//   const auto & GL = appObj.viewGhostLeft();
//   const auto & GF = appObj.viewGhostFront();
//   const auto & GR = appObj.viewGhostRight();
//   const auto & GB = appObj.viewGhostBack();

//   // auto printLam = [](std::string s,  int it, auto const & M){
//   //   std::cout << '\n';
//   //   std::cout << s << ": ";
//   //   for (int j=0; j< M.cols(); ++j){
//   //     std::cout << M(it,j) << " ";
//   //   }
//   // };

//   const auto & graph = meshObj.graph();
//   const auto & rowsBd = meshObj.graphRowsOfCellsNearBd();
//   for (decltype(rowsBd.size()) it=0; it<rowsBd.size(); ++it){
//     int smPt = rowsBd[it];
//     const auto cGID = graph(smPt, 0);
//     std::cout << cGID << ": ";

//     printLam("left   ", it, GL);
//     printLam("back   ", it, GB);
//     printLam("front  ", it, GF);
//     printLam("right  ", it, GR);

//     std::cout << "\n";
//     std::cout << "\n";
//   }

//   return true;
// }

int main()
{
  // test_sample_s7();
  // std::cout << "\n";
  // std::cout << "\n";

  if (!test_sample_s5()){
    std::cout << "\nFAILED\n";
  }
  else{
    std::cout << "PASS\n";
  }

  return 0;
}
