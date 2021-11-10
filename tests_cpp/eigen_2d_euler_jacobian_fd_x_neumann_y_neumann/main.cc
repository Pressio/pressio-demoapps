
#include "pressiodemoapps/euler2d.hpp"
#include <iomanip>
#include <array>

template<class T>
void writeToFile(const T& obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (size_t i=0; i<obj.size(); i++){
    file << std::setprecision(14) << obj(i) << " \n";
  }
  file.close();
}

template<class state_type, class mesh_t>
void modify_state(state_type & state,
		  const mesh_t & meshObj)
{

  using scalar_type = typename mesh_t::scalar_t;
  const auto gamma = 1.4;
  constexpr int numDofPerCell = 4;
  constexpr auto zero = static_cast<scalar_type>(0);
  constexpr auto one  = static_cast<scalar_type>(1);
  constexpr auto two = static_cast<scalar_type>(2);

  const auto & x= meshObj.viewX();
  const auto & y= meshObj.viewY();
  const auto gammaMinusOne = gamma - one;
  const auto gammaMinusOneInv = one/gammaMinusOne;
  constexpr auto x0 = one/two;
  constexpr auto y0 = one/two;

  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto pert = 0.001*std::sin(8.*M_PI*x(i))*std::sin(8.*M_PI*y(i));

      if (x(i) >= x0 and y(i) >= y0){
	prim = {0.5313, zero, zero, 0.4};
      }
      else if (x(i) < x0 and y(i) >= y0){
	prim = {one, 0.7276, zero, one};
      }
      else if (x(i) < x0 and y(i) < y0){
	prim = {0.8, zero, zero, one};
      }
      else if (x(i) > x0 and y(i) < y0){
	prim = {one, zero, 0.7276, one};
      }

      prim[0] += pert;
      prim[1] += pert;
      prim[2] += pert;
      prim[3] += pert;

      const auto ind = i*numDofPerCell;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = prim[0]*prim[2];
      state(ind+3) = pressiodemoapps::eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
    }
}

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
  constexpr auto order   = pda::InviscidFluxReconstruction::FirstOrder;

  // here we want to test FD Jacobian when we have Neumann BC on all boundaries
  // we can use Riemann for this since that problem has the correct BCs

  auto appObj = pda::createProblemEigen(meshObj, pda::Euler2d::Riemann, order);
  using app_t = decltype(appObj);
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  using app_jacob_t	= typename app_t::jacobian_type;

  app_state_t state = appObj.initialCondition();
  modify_state(state, meshObj);
  writeToFile(state, "IC.txt");

  auto velo = appObj.createVelocity();
  auto J = appObj.createJacobian();

  appObj.velocity(state, 0., velo);
  appObj.jacobian(state, 0., J);
  Eigen::VectorXd a = Eigen::VectorXd::Random(state.size());

  // Riemann is a non-dimensional problem, so we can use the same epsilon
  const double eps = 1e-8;

  // first order
  auto state2 = state+eps*a;
  app_rhs_t velo2(velo.size());
  appObj.velocity(state2, 0., velo2);

  auto Ja = J*a;
  auto Ja_fd = (velo2 - velo)/eps;
  for (int i=0; i<Ja.size(); ++i)
  {
    const auto diff = std::abs(Ja(i)- Ja_fd(i));

    printf(" i=%2d J(i)=%10.6f J_fd(i)=%10.6f diff=%e \n",
	   i, Ja(i), Ja_fd(i), diff);

    if (diff > 1e-4){
      std::puts("FAILED");
      return 0;
    }
  }

  // second order
  auto state3 = state-eps*a;
  app_rhs_t velo3(velo.size());
  appObj.velocity(state3, 0., velo3);
  auto Ja_fd_2 = (velo2 - velo3)/(2.*eps);
  for (int i=0; i<Ja.size(); ++i){
    const auto diff = std::abs(Ja(i)- Ja_fd_2(i));
    if (diff > 1e-6){
      std::puts("FAILED");
      return 0;
    }
  }

  std::puts("PASS");
  return 0;
}
