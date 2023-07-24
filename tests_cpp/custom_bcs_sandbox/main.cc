#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "../observer.hpp"

struct Func{
  template<class ConnecRowType, class StateT, class T>
  void operator()(const int k, /*index enumerating cells on boundary*/
		  ConnecRowType const & connectivityRow,
		  const double cellX, const double cellY,
		  const double cellWidth,
		  const StateT & currentState,
		  int numDofPerCell,
		  pressiodemoapps::GhostRelativeLocation ghostRelLoc,
		  T & ghostValues) const
  {
    /* here we use kind of arbitrary numbers just to prove it works*/

    assert(ghostValues.size() == numDofPerCell);

    const int cellGID = connectivityRow[0];
    const auto uIndex  = cellGID*numDofPerCell;

    // Dirichlet for dof0
    ghostValues[0] = 0.2;

    // non-homog Neumann for dof1
    ghostValues[1] = cellWidth*0.01 + currentState(uIndex+1);

    // homog Neumann for dof2,3
    ghostValues[2] = currentState(uIndex+2);
    ghostValues[3] = currentState(uIndex+3);
  }

  template<class ConnecRowType, class FactorsType>
  void operator()(ConnecRowType const & connectivityRow,
		  const double cellX, const double cellY,
		  int numDofPerCell,
		  int axis, /* =1 for x, =2 for y*/
		  FactorsType & factorsForBCJac) const
  {
    factorsForBCJac = {0.,        // 0 should always be set for Dirichlet
		       1.,1.,1.}; // 1 should always be used for Neumann
  }
};

template<class AppT>
bool verifyJacobian(AppT & appObj)
{
  auto state = appObj.initialCondition();
  auto velo = appObj.createRightHandSide();
  auto J = appObj.createJacobian();

  appObj.rightHandSide(state, 0., velo);
  appObj.jacobian(state, 0., J);
  Eigen::VectorXd a = Eigen::VectorXd::Random(state.size());
  auto Ja = J*a;

  // Riemann is a non-dimensional problem, so we can use the same epsilon
  const double eps = 1e-8;

  auto verify = [&](const auto & Ja_fd, double tol){
    for (int i=0; i< Ja.size(); ++i){
      const auto diff = std::abs(Ja(i)- Ja_fd(i));
      printf(" i=%2d J(i)=%10.6f J_fd(i)=%10.6f diff=%e \n",
	     i, Ja(i), Ja_fd(i), diff);
      if (diff > tol || Ja(i)==0.){ return false; }
    }
    return true;
  };

  auto state2 = state+eps*a;
  decltype(velo) velo2(velo.size());
  appObj.rightHandSide(state2, 0., velo2);
  auto Ja_fd = (velo2 - velo)/eps;
  bool res1  = verify(Ja_fd, 1e-5);


  // second order FD
  auto state3 = state-eps*a;
  decltype(velo) velo3(velo.size());
  appObj.rightHandSide(state3, 0., velo3);
  auto Ja_fd2 = (velo2 - velo3)/(2.*eps);
  bool res2 = verify(Ja_fd2, 5e-7);

  return res1 && res2;
}

int main()
{
  namespace pda = pressiodemoapps;

  const auto inviscidScheme   = pda::InviscidFluxReconstruction::FirstOrder;
  const auto probId  = pda::Euler2d::RiemannCustomBCs;
  const int initCond = 2;

  {
    const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("./fullmesh");
    Func bcF;
    auto appObj = pda::create_problem_eigen(meshObj, probId, inviscidScheme, bcF, initCond);
    bool res = verifyJacobian(appObj);
    if (!res){
      std::puts("\nFAILED");
      return 0;
    }
  }

  {
    const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("./samplemesh");
    Func bcF;
    auto appObj = pda::create_problem_eigen(meshObj, probId, inviscidScheme, bcF, initCond);
    bool res = verifyJacobian(appObj);
    if (!res){
      std::puts("\nFAILED");
      return 0;
    }
  }

  std::puts("PASS");
  return 0;
}
