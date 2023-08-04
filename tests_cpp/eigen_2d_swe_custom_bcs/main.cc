#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/swe2d.hpp"
#include "../observer.hpp"

struct Dirichlet
{
  template<class ConnecRowType, class StateT, class T>
  void operator()(const int /*unused*/, ConnecRowType const & connectivityRow,
		  const double cellX, const double cellY,
		  const StateT & currentState, int numDofPerCell,
		  const double cellWidth, T & ghostValues) const
  {
    assert(ghostValues.size() == numDofPerCell);
    assert(numDofPerCell == 3);

    ghostValues(0) = 0.00001;
    ghostValues(1) = 0.004;
    ghostValues(2) = 0.001;
  }

  template<class ConnecRowType, class FactorsType>
  void operator()(ConnecRowType const & connectivityRow,
		  const double cellX, const double cellY,
		  int numDofPerCell, FactorsType & factorsForBCJac) const
  {
    assert(numDofPerCell == 3);
    factorsForBCJac = {0.,0.,0.};
  }
};

struct HomogNeumann
{
  template<class ConnecRowType, class StateT, class T>
  void operator()(const int /*unused*/, ConnecRowType const & connectivityRow,
		  const double cellX, const double cellY,
		  const StateT & currentState, int numDofPerCell,
		  const double cellWidth, T & ghostValues) const
  {
    assert(ghostValues.size() == numDofPerCell);
    assert(numDofPerCell == 3);

    const int cellGID = connectivityRow[0];
    const auto uIndex  = cellGID*numDofPerCell;
    ghostValues[0] = currentState(uIndex);
    ghostValues[1] = currentState(uIndex+1);
    ghostValues[2] = currentState(uIndex+2);
  }

  template<class ConnecRowType, class FactorsType>
  void operator()(ConnecRowType const & connectivityRow,
		  const double cellX, const double cellY,
		  int numDofPerCell, FactorsType & factorsForBCJac) const
  {
    assert(numDofPerCell == 3);
    factorsForBCJac = {1.,1.,1.};
  }
};

template<class AppT, class MeshT>
bool verifyJacobian(AppT & appObj, MeshT & meshObj)
{
  auto state = appObj.initialCondition();
  // perturb initial condition
  for (int i=0; i<state.size(); ++i){
    state(i) += 0.5;
  }
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

  std::cout << " doing first order FD \n" ;
  auto state2 = state+eps*a;
  decltype(velo) velo2(velo.size());
  appObj.rightHandSide(state2, 0., velo2);
  auto Ja_fd = (velo2 - velo)/eps;
  bool res1  = verify(Ja_fd, 1e-5);

  // second order FD
  std::cout << " doing second order FD \n" ;
  auto state3 = state-eps*a;
  decltype(velo) velo3(velo.size());
  appObj.rightHandSide(state3, 0., velo3);
  auto Ja_fd2 = (velo2 - velo3)/(2.*eps);
  bool res2 = verify(Ja_fd2, 1e-6);

  std::cout << "\n verifyJacobian FINISHED \n";
  return res1 && res2;
}

int main()
{
  namespace pda = pressiodemoapps;

  const auto inviscidScheme   = pda::InviscidFluxReconstruction::FirstOrder;
  const auto probId  = pda::Swe2d::CustomBCs;

  {
    std::cout << " run on full mesh \n" ;
    const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("./fullmesh");
    auto appObj = pda::create_problem_eigen(meshObj, probId, inviscidScheme,
					    Dirichlet(), Dirichlet(),
					    HomogNeumann(), HomogNeumann());

    bool res = verifyJacobian(appObj, meshObj);
    if (!res){
      std::puts("\nFAILED");
      return 0;
    }
  }

  {
    std::cout << " \nrun on sample mesh \n" ;
    const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("./samplemesh");
    auto appObj = pda::create_problem_eigen(meshObj, probId, inviscidScheme,
					    Dirichlet(), Dirichlet(),
					    HomogNeumann(), HomogNeumann());
    bool res = verifyJacobian(appObj, meshObj);
    if (!res){
      std::puts("\nFAILED");
      return 0;
    }
  }

  std::puts("PASS");
  return 0;
}
