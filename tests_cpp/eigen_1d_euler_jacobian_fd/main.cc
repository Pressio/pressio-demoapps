
#include "pressio_containers.hpp"
#include "euler1d.hpp"
#include <iomanip>

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
  constexpr int numDofPerCell = 3;
  constexpr auto zero = static_cast<scalar_type>(0);
  constexpr auto one  = static_cast<scalar_type>(1);

  const auto & x= meshObj.viewX();
  std::array<scalar_type, numDofPerCell> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto pert = 0.01*std::sin(8.*M_PI*x(i));

      if (x(i) <= zero){
	prim[0] = one + pert;
	prim[1] = zero + pert;
	prim[2] = one + pert;
      }

      if (x(i) > zero){
	prim[0] = static_cast<scalar_type>(0.125) + pert;
	prim[1] = zero + pert;
	prim[2] = static_cast<scalar_type>(0.1) + pert;
      }

      const auto ind = i*numDofPerCell;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = pressiodemoapps::eulerEquationsComputeEnergyFromPrimitive(gamma, prim);
    }
}

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
  constexpr auto order   = pda::InviscidFluxReconstruction::FirstOrder;

  // even if we set Sod, this does not really matter
  auto appObj = pda::createProblemEigen(meshObj, pda::Euler1d::Sod, order);
  using app_t = decltype(appObj);
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  using app_jacob_t	= typename app_t::jacobian_type;

  using ode_state_t = pressio::containers::Vector<app_state_t>;
  using ode_res_t   = pressio::containers::Vector<app_rhs_t>;
  using ode_jac_t   = pressio::containers::SparseMatrix<app_jacob_t>;

  app_state_t state = appObj.initialCondition();
  modify_state(state, meshObj);
  writeToFile(state, "IC.txt");

  auto velo = appObj.createVelocity();
  auto J = appObj.createJacobian();

  const double eps = 1e-8;

  // make sure repeated evaluations work
  // not just a single time
  for (int i=0; i<5; ++i)
  {
    appObj.velocity(state, 0., velo);
    appObj.jacobian(state, 0., J);

    Eigen::VectorXd a = Eigen::VectorXd::Random(state.size());

    // first order
    auto state2 = state+eps*a;
    app_rhs_t velo2(velo.size());
    appObj.velocity(state2, 0., velo2);

    auto Ja = J*a;
    auto Ja_fd = (velo2 - velo)/eps;
    for (int i=0; i<Ja.size(); ++i)
    {
      const auto diff = std::abs(Ja(i)- Ja_fd(i));
      std::cout << i << " "
		<< diff << " "
		<< Ja(i) << " "
		<< Ja_fd(i) << " "
		<< std::endl;

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
  }

  std::puts("PASS");
  return 0;
}
