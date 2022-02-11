
#include "pressiodemoapps/euler2d.hpp"
#include <iomanip>

template<class state_type, class mesh_t>
void modify_state(state_type & state, const mesh_t & meshObj)
{
  const auto & x= meshObj.viewX();
  const auto & y= meshObj.viewY();
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto pert = 0.01*std::sin(24.*M_PI*x(i)*y(i));
      const auto ind = i*4;
      state[ind] += pert;
      state[ind+1] += pert;
      state[ind+2] += pert;
      state[ind+3] += pert;
    }
}

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");

  const auto scheme = ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder;
  const auto probId = pda::Euler2d::RayleighTaylor;
  const double amplitude = 0.025;
  auto appObj     = pda::create_problem_eigen(meshObj, probId, scheme, amplitude);
  using app_t = decltype(appObj);
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  using app_jacob_t	= typename app_t::jacobian_type;

  app_state_t state = appObj.initialCondition();
  modify_state(state, meshObj);

  auto velo = appObj.createVelocity();
  auto J = appObj.createJacobian();

  const double eps = 1e-8;

  // make sure repeated evaluations work
  // not just a single time
  for (int loop=0; loop<5; ++loop)
  {
    appObj.jacobian(state, 0., J);
    appObj.velocity(state, 0., velo);

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

      if (loop==0){
        printf(" i=%2d J(i)=%10.6f J_fd(i)=%10.6f diff=%e \n",
	       i, Ja(i), Ja_fd(i), diff);
      }

      if (diff > 5e-5){
        std::puts("FAILED");
        return 0;
      }
    }
  }

  std::puts("PASS");
  return 0;
}
