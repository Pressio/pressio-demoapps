
#include "pressiodemoapps/advection_diffusion2d.hpp"
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

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");

  const auto ischeme  = pda::InviscidFluxReconstruction::Weno5;
  const auto vscheme  = pda::ViscousFluxReconstruction::FirstOrder;
  const auto probId  = pda::AdvectionDiffusion2d::BurgersPeriodic;
  auto appObj        = pda::create_problem_eigen(meshObj, probId, ischeme, vscheme);
  using app_t = decltype(appObj);
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  using app_jacob_t	= typename app_t::jacobian_type;

  app_state_t state = appObj.initialCondition();
  writeToFile(state, "IC.txt");

  auto velo = appObj.createVelocity();
  auto J = appObj.createJacobian();

  const double eps = 1e-8;

  // make sure repeated evaluations work
  for (int loop=0; loop<6; ++loop)
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
      if (loop==0){
	printf(" i=%2d J(i)=%10.6f J_fd(i)=%10.6f diff=%e \n",
	       i, Ja(i), Ja_fd(i), diff);
      }

      if (diff > 5e-5){
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

      if (loop==0){
	printf(" i=%2d J(i)=%10.6f J_fd(i)=%10.6f diff=%e \n",
	       i, Ja(i), Ja_fd(i), diff);
      }

      if (diff > 1e-6){
	std::puts("FAILED");
	return 0;
      }
    }

  }

  std::puts("PASS");
  return 0;
}
