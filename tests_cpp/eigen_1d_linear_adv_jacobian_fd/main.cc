
#include "pressio/type_traits.hpp"
#include "pressiodemoapps/advection.hpp"
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
  const auto & x = meshObj.viewX();
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto pert = 0.1*std::sin(16.*M_PI*x(i));
      state(i) += pert;
    }
}

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");

  constexpr auto order   = pda::InviscidFluxReconstruction::FirstOrder;
  const auto probid = pda::Advection1d::PeriodicLinear;
  auto appObj       = pda::createProblemEigen(meshObj, probid, order);

  using app_t = decltype(appObj);
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  using app_jacob_t	= typename app_t::jacobian_type;

  app_state_t state = appObj.initialCondition();
  modify_state(state, meshObj);
  writeToFile(state, "IC.txt");

  const double eps = 1e-8;
  auto velo = appObj.createVelocity();
  auto J = appObj.createJacobian();

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
      printf(" i=%2d J(i)=%10.6f J_fd(i)=%10.6f diff=%e \n", i, Ja(i), Ja_fd(i), diff);

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
