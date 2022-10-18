
#include "pressio/type_traits.hpp"
#include "pressiodemoapps/advection1d.hpp"
#include <iomanip>
#include <random>

template<class T>
void writeToFile(const T& obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (decltype(obj.size()) i=0; i<obj.size(); i++){
    file << std::setprecision(14) << obj(i) << " \n";
  }
  file.close();
}

template<class state_type, class mesh_t>
void modify_state(state_type & state,
		  const mesh_t & meshObj)
{
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<> dist(0, 1);

  const auto & x = meshObj.viewX();
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      //const auto pert = 2.*std::sin(16.*M_PI*(x(i)-0.2));
      state(i) = dist(rng);//pert;
    }
}

int main()
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");

#ifdef USE_WENO5
  const auto order   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#elif defined USE_FIRSTORDER
  constexpr auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  auto appObj = pda::create_linear_advection_1d_problem_eigen(
    meshObj, order, pda::InviscidFluxScheme::Rusanov, 2.0);

  using app_t = decltype(appObj);
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;

  app_state_t state = appObj.initialCondition();
  modify_state(state, meshObj);
  writeToFile(state, "IC.txt");

  const double eps = 1e-8;
  auto velo = appObj.createRightHandSide();
  auto J = appObj.createJacobian();

  // make sure repeated evaluations work
  // not just a single time
  for (int loop=0; loop<5; ++loop)
  {
    std::cout << "!!!! VELOCITY !!!\n";
    appObj.rightHandSide(state, 0., velo);
    std::cout << "!!!! JACOBIAN !!!\n";
    appObj.jacobian(state, 0., J);

    Eigen::VectorXd a = Eigen::VectorXd::Random(state.size());
    auto Ja = J*a;

    // first order
    auto state2 = state+eps*a;
    app_rhs_t velo2(velo.size());
    appObj.rightHandSide(state2, 0., velo2);

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
    appObj.rightHandSide(state3, 0., velo3);
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
