
#include "pressiodemoapps/diffusion_reaction2d.hpp"
#include <iomanip>

template<class T>
void writeToFile(const T& obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (int i=0; i<obj.size(); i++){
    file << std::setprecision(14) << obj(i) << " \n";
  }
  file.close();
}

template<class state_type, class mesh_t>
void modify_state(state_type & state,
		  const mesh_t & meshObj)
{

  const auto & x= meshObj.viewX();
  const auto & y= meshObj.viewY();
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto ind = i*2;
      const auto pert = 0.01*std::sin(24.*M_PI*x(i)*y(i));
      state[ind] += pert;
      state[ind+1] += pert;
    }
}

int main()
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");

  const auto scheme  = pda::ViscousFluxReconstruction::FirstOrder;
  const auto probId  = pda::DiffusionReaction2d::GrayScott;
  auto appObj        = pda::create_problem_eigen(meshObj, probId, scheme); //, 0.0002, 0.0002/4., 0.03, 0.062);
  using app_t = decltype(appObj);
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;

  app_state_t state = appObj.initialCondition();
  modify_state(state, meshObj);
  writeToFile(state, "IC.txt");

  auto velo = appObj.createRightHandSide();
  auto J = appObj.createJacobian();

  const double eps = 1e-8;

  // make sure repeated evaluations work
  for (int loop=0; loop<2; ++loop)
  {
    appObj.rightHandSide(state, 0., velo);
    appObj.jacobian(state, 0., J);

    Eigen::VectorXd a = Eigen::VectorXd::Random(state.size());

    // first order
    auto state2 = state+eps*a;
    app_rhs_t velo2(velo.size());
    appObj.rightHandSide(state2, 0., velo2);

    auto Ja = J*a;
    auto Ja_fd = (velo2 - velo)/eps;
    for (int i=0; i<Ja.size(); ++i)
    {
      const auto diff = std::abs(Ja(i)- Ja_fd(i));
      if (loop==0){
	       printf(" i=%2d J(i)=%10.6f J_fd(i)=%10.6f diff=%e \n",
	       i, Ja(i), Ja_fd(i), diff);
      }

      if (diff > 1e-6){
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

      if (loop==0){
	     printf(" i=%2d J(i)=%10.6f J_fd(i)=%10.6f diff=%e \n",
	       i, Ja(i), Ja_fd(i), diff);
      }

      if (diff > 1e-7){
	      std::puts("FAILED");
	      return 0;
      }
    }

  }

  std::puts("PASS");
  return 0;
}
