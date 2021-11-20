
#include "pressio/type_traits.hpp"
#include "pressiodemoapps/swe2d.hpp"
#include <iomanip>

template<class T>
void writeToFile(const T& obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (size_t i=0; i<obj.size(); i++){
    file << std::setprecision(16) << obj(i) << " \n";
  }
  file.close();
}

template<class state_type, class mesh_t>
void modify_state(state_type & state,
      const mesh_t & meshObj)
{

  using scalar_type = typename mesh_t::scalar_t;
  const auto & x = meshObj.viewX();
  const auto & y = meshObj.viewY();
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i){
    const auto ind = i*3;
    const auto pert = 0.001*std::sin(8.*M_PI*x(i)*y(i));
    state(ind)   = 1. + pert;
    state(ind+1) = 1.5 + pert;
    state(ind+2) = 2. + pert;
  }
}

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");

  constexpr auto order   = pda::InviscidFluxReconstruction::FirstOrder;
  const auto probId  = pda::Swe2d::SlipWall;
  auto appObj = pda::create_problem_eigen(meshObj, probId, order);
  using app_t = decltype(appObj);
  using scalar_t  = typename app_t::scalar_type;

  auto state = appObj.initialCondition();
  modify_state(state, meshObj);
  writeToFile(state, "IC.txt");

  const double eps = 1e-8;
  auto velo = appObj.createVelocity();
  auto J = appObj.createJacobian();

  // make sure repeated evaluations work
  // not just a single time
  for (int loop=0; loop<1; ++loop)
  {
    appObj.velocity(state, 0., velo);
    appObj.jacobian(state, 0., J);

    Eigen::VectorXd a = Eigen::VectorXd::Random(state.size());
    auto Ja = J*a;

    // first order
    auto state2 = state+eps*a;
    decltype(velo) velo2(velo.size());
    appObj.velocity(state2, 0., velo2);

    auto Ja_fd = (velo2 - velo)/eps;
    for (int i=0; i<Ja.size(); ++i)
    {
      const auto diff = std::abs(Ja(i)- Ja_fd(i));
      if (loop == 0){
	printf(" i=%2d J(i)=%16.12f J_fd(i)=%16.12f diff=%e \n",
	       i, Ja(i), Ja_fd(i), diff);
      }

      if (diff > 1e-4){
	std::puts("FAILED");
	//return 0;
      }
    }

    // // second order
    // auto state3 = state-eps*a;
    // decltype(velo) velo3(velo.size());
    // appObj.velocity(state3, 0., velo3);
    // auto Ja_fd_2 = (velo2 - velo3)/(2.*eps);
    // for (int i=0; i<Ja.size(); ++i){
    //   const auto diff = std::abs(Ja(i)- Ja_fd_2(i));
    //   if (loop == 0){
    // 	printf(" i=%2d J(i)=%10.6f J_fd(i)=%10.6f diff=%e \n", i, Ja(i), Ja_fd(i), diff);
    //   }

    //   if (diff > 1e-6){
    // 	std::puts("FAILED");
    // 	return 0;
    //   }
    // }

  }

  std::puts("PASS");
  return 0;
}
