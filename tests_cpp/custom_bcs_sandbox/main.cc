#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "../observer.hpp"

struct Func{

  /* this overload is called to set the ghost values */
  template<class ConnecRowType, class StateT, class T>
  void operator()(const int k,
		  /*
		    k:
		    index enumerating the cells on boundary: 0, 1, 2, .., N
		    where N is the number of cells that are near the boundary
		    as per the definition of the mesh object.
		    Not sure if we need this but just in case.
		  */
		  ConnecRowType const & connectivityRow,
		  /*connectivityRow:
		    contains the row of the connectivity graph of the cell being handled */

		  const double cellX, const double cellY,
		  /* cellX, cellY: coordinates of the cell being handled */

		  const StateT & currentState,
		  /* currentState:
		     a const ref to the state vector currently used inside the app object */

		  int numDofPerCell,
		  /*  numDofPerCell: self-explanatory */

		  int workingAxis,
		  /* can be 1 or 2:
		     - if workingAxis=1, we are dealing with ghosts alogn x
		     - if workingAxis=2, dealing with ghosts along y*/

		  pressiodemoapps::GhostRelativeLocation ghostRelLoc,
		  /* ghostRelLoc:
		     this tells you the position of the ghost relative to the boundary cell we are dealing,
		     it can have one of the values: Left,Right,Back,Front */

		  const double cellWidth,
		  /* cellWidth:
		     - if workingAxis==1, then cellWidth = dx
		     - if workingAxis==2, then cellWidth = dy */

		  T & ghostValues
		  /* the "vector" of ghost values to overwrite for all dofs */

		  ) const
  {

    // preconditions
    assert(ghostValues.size() == numDofPerCell);
    if (workingAxis==1){
      assert(ghostRelLoc == pressiodemoapps::GhostRelativeLocation::Left
	     || ghostRelLoc == pressiodemoapps::GhostRelativeLocation::Right);
    }
    if (workingAxis==2){
      assert(ghostRelLoc == pressiodemoapps::GhostRelativeLocation::Back
	     || ghostRelLoc == pressiodemoapps::GhostRelativeLocation::Front);
    }

    const int cellGID = connectivityRow[0];
    const auto uIndex  = cellGID*numDofPerCell;

    /*
      IMPORTANT: this code below is totally random,
      we use kind of arbitrary numbers just to prove it works.
      Just for testing, on left, back and right boundary we do something,
      while on the top boundary we do something else.

      Also, we set different BCs for each dofs, and we know here we
      have 4 dofs because we are doing 2d Euler
    */

    if (ghostRelLoc == pressiodemoapps::GhostRelativeLocation::Left){
      ghostValues[0] = 0.2; // Dirichlet
      ghostValues[1] = cellWidth*0.01 + currentState(uIndex+1); // non-homog Neumann
      ghostValues[2] = currentState(uIndex+2); //homog Neumann
      ghostValues[3] = currentState(uIndex+3); //homog Neumann
    }

    else if (ghostRelLoc == pressiodemoapps::GhostRelativeLocation::Back){
      ghostValues[0] = cellWidth*0.01 + currentState(uIndex); // non-homog Neumann
      ghostValues[1] = 0.0001; // Dirichlet
      ghostValues[2] = currentState(uIndex+2); //homog Neumann
      ghostValues[3] = currentState(uIndex+3); //homog Neumann
    }

    else{
      for (int i=0; i<numDofPerCell; ++i){
	ghostValues[i] = currentState(uIndex+i); // homog Neumann for all
      }
    }
  }

  template<class ConnecRowType, class FactorsType>
  void operator()(ConnecRowType const & connectivityRow,
		  const double cellX, const double cellY,
		  int numDofPerCell,
		  int workingAxis,
		  pressiodemoapps::GhostRelativeLocation ghostRelLoc,
		  FactorsType & factorsForBCJac) const
  {
    // 0 should always be set for Dirichlet
    // 1 should always be used for Neumann

    if (ghostRelLoc == pressiodemoapps::GhostRelativeLocation::Left){
      factorsForBCJac = {0.,1.,1.,1.};
    }
    else if (ghostRelLoc == pressiodemoapps::GhostRelativeLocation::Back){
      factorsForBCJac = {1.,0.,1.,1.};
    }
    else{
      factorsForBCJac = {1.,1.,1.,1.};
    }
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
  bool res2 = verify(Ja_fd2, 5e-7);

  std::cout << "\n verifyJacobian FINISHED \n";
  return res1 && res2;
}

int main()
{
  namespace pda = pressiodemoapps;

  const auto inviscidScheme   = pda::InviscidFluxReconstruction::FirstOrder;
  const auto probId  = pda::Euler2d::RiemannCustomBCs;
  const int initCond = 2;

  {
    std::cout << " run on full mesh \n" ;
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
    std::cout << " \nrun on sample mesh \n" ;
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
