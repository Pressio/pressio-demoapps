
#include "pressiodemoapps/advection_diffusion_reaction2d.hpp"
#include <iomanip>

template<class T>
void writeToFile(const T& obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (int i=0; i<obj.size(); i++){
    file << std::setprecision(14) << obj(i) << " \n";
  }
  file.close();
}

template<class state_type, class mesh_t>
void modify_state(state_type & state, const mesh_t & meshObj){
  const auto & x= meshObj.viewX();
  const auto & y= meshObj.viewY();
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i){
      const auto ind = i;
      const auto pert = 0.01*std::sin(24.*M_PI*x(i)*y(i));
      state[ind] += pert;
    }
}

template<class MeshObjT, class AppObjT, class T>
void runWithSpecialTreatmentAtBd(const MeshObjT & meshObj,
				 AppObjT & appObjFirstOrder,
				 AppObjT & appObj,
				 const typename AppObjT::state_type & state,
				 T eps,
				 T tol1FD,
				 T tol2FD,
				 int numLoops)
{
  // for near the boundaries, we use first-order jacobians

  using app_rhs_t = typename AppObjT::velocity_type;

  auto J     = appObj.createJacobian();
  auto velo  = appObj.createRightHandSide();
  auto velo2 = appObj.createRightHandSide();
  auto velo3 = appObj.createRightHandSide();

  for (int loop=0; loop<numLoops; ++loop)
  {
    Eigen::VectorXd a = Eigen::VectorXd::Random(state.size());

    // compute Jacobian for the appObj of interest
    appObj.rightHandSide(state, 0., velo);
    appObj.jacobian(state, 0., J);
    auto Ja = J*a;

    // prepare for FD
    auto statePlusEps  = state+eps*a;
    auto stateMinusEps = state-eps*a;

    auto verify = [&](auto cells, const app_rhs_t & Ja_fd_in, T tol)
    {
      for (decltype(cells.size()) i=0; i<cells.size(); ++i){
	const int row = cells[i];
	const auto diff = std::abs(Ja(row) - Ja_fd_in(row));
	  if (loop==0){
	    printf(" row=%2d J(row)=%10.6f J_fd(row)=%10.6f diff=%e \n",
		   row, Ja(row), Ja_fd_in(row), diff);
	  }
	  if (diff > tol || Ja_fd_in(row)==0.){ return false; }
      }
      return true;
    };

    //
    // 1.
    // verify that near the boundaries the Jacobian
    // matches the first-order one

    // compute J for first-order app obj
    auto velo_fo  = appObj.createRightHandSide();
    auto velo_fo2 = appObj.createRightHandSide();
    appObjFirstOrder.rightHandSide(state, 0., velo_fo);
    appObjFirstOrder.rightHandSide(statePlusEps, 0., velo_fo2);
    auto Ja_fd_fo = (velo_fo2 - velo_fo)/eps;
    const auto & graphRowsNearBd = meshObj.graphRowsOfCellsNearBd();
    bool res = verify(graphRowsNearBd, Ja_fd_fo, tol1FD);
    if (!res){ std::puts("FAILED"); }

    std::cout << "\n";
    std::cout << "\n";

    auto velo_fo3 = appObj.createRightHandSide();
    appObjFirstOrder.rightHandSide(stateMinusEps, 0., velo_fo3);
    auto Ja_fd_fo2 = (velo_fo2 - velo_fo3)/(2.*eps);
    res = verify(graphRowsNearBd, Ja_fd_fo2, tol2FD);
    if (!res){ std::puts("FAILED"); }

    std::cout << "\n";
    std::cout << "\n";

    //
    // 2.
    // verify that for internal cells the Jacobian matches the FD one
    appObj.rightHandSide(state, 0., velo);
    appObj.rightHandSide(statePlusEps, 0., velo2);
    auto Ja_fd = (velo2 - velo)/eps;
    const auto & graphRowsInternal = meshObj.graphRowsOfCellsAwayFromBd();
    bool res2 = verify(graphRowsInternal, Ja_fd, tol1FD);
    if (!res2){ std::puts("FAILED"); }

    std::cout << "\n";
    std::cout << "\n";

    auto velo3 = appObj.createRightHandSide();
    appObj.rightHandSide(stateMinusEps, 0., velo3);
    auto Ja_fd2 = (velo2 - velo3)/(2.*eps);
    res = verify(graphRowsInternal, Ja_fd2, tol2FD);
    if (!res){ std::puts("FAILED"); }
  }
}

template<class AppObjT, class T>
void defaultCase(AppObjT & appObj,
		const typename AppObjT::state_type & state,
		T eps,
		T tol1FD,
		T tol2FD,
		int numLoops)
{
  using app_rhs_t = typename AppObjT::velocity_type;

  auto J = appObj.createJacobian();
  auto velo = appObj.createRightHandSide();
  auto velo2 = appObj.createRightHandSide();

  for (int loop=0; loop<numLoops; ++loop)
  {
    Eigen::VectorXd a = Eigen::VectorXd::Random(state.size());

    appObj.jacobian(state, 0., J);
    auto Ja = J*a;

    auto verify = [&](const app_rhs_t & Ja_fd_in, T tol) -> bool
    {
      for (int i=0; i<Ja.size(); ++i){
	const auto diff = std::abs(Ja(i)- Ja_fd_in(i));
	if (loop==0){
	  printf(" i=%2d J(i)=%10.6f J_fd(i)=%10.6f diff=%e \n",
		 i, Ja(i), Ja_fd_in(i), diff);
	}

	if (diff > tol || Ja_fd_in(i)==0.){ return false; }
      }
      return true;
    };

    //
    // first order FD
    //
    appObj.rightHandSide(state, 0., velo);
    auto state2 = state+eps*a;
    appObj.rightHandSide(state2, 0., velo2);
    auto Ja_fd = (velo2 - velo)/eps;
    bool res = verify(Ja_fd, tol1FD);
    if (!res){ std::puts("FAILED"); }

    std::cout << "\n";
    std::cout << "\n";

    //
    // second order FD
    //
    auto state3 = state-eps*a;
    app_rhs_t velo3(velo.size());
    appObj.rightHandSide(state3, 0., velo3);
    auto Ja_fd_2 = (velo2 - velo3)/(2.*eps);
    bool res2 = verify(Ja_fd_2, tol2FD);
    if (!res2){ std::puts("FAILED"); }
  }
}


int main()
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");

#ifdef USE_WENO5
  const auto inviscidScheme = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  const auto inviscidScheme = pda::InviscidFluxReconstruction::Weno3;
#elif defined USE_FIRSTORDER
  const auto inviscidScheme = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probId  = pda::AdvectionDiffusionReaction2d::ProblemA;
  auto appObj        = pda::create_problem_eigen(meshObj, probId, inviscidScheme);

  auto state = appObj.initialCondition();
  modify_state(state, meshObj);
  writeToFile(state, "IC.txt");

  constexpr int numLoops = 1;
  constexpr double eps = 1e-8;
#if defined USE_FIRSTORDER
  constexpr double tolForFirstOrderFD = 1e-7;
  constexpr double tolForSecondOrderFD = 5e-8;
#else
  constexpr double tolForFirstOrderFD = 5e-4;
  constexpr double tolForSecondOrderFD = 5e-8;
#endif

  /* need to make a distinction because we know that
     - when inviscid is first order, everything is first order
     - when inviscid is weno3 or 5, near the boundaries we still use first order jacobians
   */
  if (inviscidScheme == pda::InviscidFluxReconstruction::FirstOrder){
    defaultCase(appObj, state, eps, tolForFirstOrderFD,
		tolForSecondOrderFD, numLoops);
  }

  else{
    const auto inviscidSchemeFo = pda::InviscidFluxReconstruction::FirstOrder;
    auto appObjFirstOrder = pda::create_problem_eigen(meshObj, probId,  inviscidSchemeFo);
    runWithSpecialTreatmentAtBd(meshObj, appObjFirstOrder, appObj, state, eps,
				tolForFirstOrderFD, tolForSecondOrderFD, numLoops);
  }

  std::puts("PASS");
  return 0;
}
