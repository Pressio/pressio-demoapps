#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "../observer.hpp"

struct Dirichlet
{
  mutable std::array<double, 4> prim = {};
  mutable std::array<double, 4> cons = {};

  template<class ConnecRowType, class StateT, class T>
  void operator()(const int /*unused*/, ConnecRowType const & connectivityRow,
		  const double cellX, const double cellY,
		  const StateT & currentState, int numDofPerCell,
		  const double cellWidth, T & ghostValues) const
  {
    assert(ghostValues.size() == numDofPerCell);
    const auto x0 = static_cast<double>(0.8);
    const double gammaMinusOne = 1.4 - 1;
    const double gammaMinusOneInv = 1/gammaMinusOne;

    if (cellX >= x0){
      prim = {1.5, 0., 0., 1.5};
    }
    else{
      prim = {0.5323, 1.206, 0., 0.3};
    }

    ghostValues(0) = prim[0];
    ghostValues(1) = prim[0]*prim[1];
    ghostValues(2) = prim[0]*prim[2];
    ghostValues(3) = pressiodemoapps::eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
  }

  template<class ConnecRowType, class FactorsType>
  void operator()(ConnecRowType const & connectivityRow,
		  const double cellX, const double cellY,
		  int numDofPerCell, FactorsType & factorsForBCJac) const
  {
    factorsForBCJac = {0.,0.,0.,0.};
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
    const int cellGID = connectivityRow[0];
    const auto uIndex  = cellGID*numDofPerCell;
    ghostValues[0] = currentState(uIndex);
    ghostValues[1] = currentState(uIndex+1);
    ghostValues[2] = currentState(uIndex+2);
    ghostValues[3] = currentState(uIndex+3);
  }

  template<class ConnecRowType, class FactorsType>
  void operator()(ConnecRowType const & connectivityRow,
		  const double cellX, const double cellY,
		  int numDofPerCell, FactorsType & factorsForBCJac) const
  {
    factorsForBCJac = {1.,1.,1.,1.};
  }
};

template<class state_type, class mesh_t>
void my_initial_condition(state_type & state,
			  const mesh_t & meshObj)
{
  using scalar_type = double;

  constexpr int numDofPerCell = 4;
  constexpr auto zero = static_cast<scalar_type>(0);
  constexpr auto one = static_cast<scalar_type>(1);
  const auto x0 = static_cast<scalar_type>(0.5);
  const auto gammaMinusOne = 1.4 - one;
  const auto gammaMinusOneInv = one/gammaMinusOne;

  const auto & x= meshObj.viewX();
  std::array<scalar_type, numDofPerCell> prim = {};
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      if (x(i) >= x0){
	prim = {1.5, zero, zero, 1.5};
      }
      else{
	prim = {0.5323, 1.206, zero, 0.3};
      }

      const auto ind = i*numDofPerCell;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = prim[0]*prim[2];
      state(ind+3) = pressiodemoapps::eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
    }
}

template<class AppT, class MeshT>
bool verifyJacobian(AppT & appObj, MeshT & meshObj)
{
  auto state = appObj.initialCondition();
  my_initial_condition(state, meshObj);

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
  bool res1  = verify(Ja_fd, 5e-6);

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
  const auto probId  = pda::Euler2d::RiemannCustomBCs;
  const int initCond = 2;

  Dirichlet bcF1;
  HomogNeumann bcF2;

  {
    std::cout << " run on full mesh \n" ;
    const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("./fullmesh");
    auto appObj = pda::create_problem_eigen(meshObj, probId, inviscidScheme,
					    bcF1, bcF1, bcF2, bcF2, initCond);

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
					    bcF1, bcF1, bcF2, bcF2, initCond);
    bool res = verifyJacobian(appObj, meshObj);
    if (!res){
      std::puts("\nFAILED");
      return 0;
    }
  }

  std::puts("PASS");
  return 0;
}
