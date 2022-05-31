
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/advection_diffusion2d.hpp"
#include "../observer.hpp"

template<class T>
void writeToFileRank1(const T & obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (size_t i=0; i<obj.rows(); i++){
    file << std::setprecision(14) << obj(i) << " \n";
  }
  file.close();
}

template<class T>
void writeToFileSparseMat(const T & obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (size_t i=0; i<obj.rows(); i++){
    for (size_t j=0; j<obj.cols(); j++){
      file << std::setprecision(14) << obj.coeff(i,j) << " ";
    }
    file << " \n";
  }
  file.close();
}

template<class state_type, class mesh_t>
void modify_state(state_type & state,
      const mesh_t & meshObj)
{

  using scalar_type = typename mesh_t::scalar_t;
  constexpr auto zero = static_cast<scalar_type>(0);
  constexpr auto one  = static_cast<scalar_type>(1);

  const auto & x= meshObj.viewX();
  const auto & y= meshObj.viewY();
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto ind = i*2;
      const auto pert = 0.1*std::sin(24.*M_PI*x(i)*y(i));
      state[ind] += pert;
      state[ind+1] += pert;
    }
}

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
  const auto probid = pda::AdvectionDiffusion2d::BurgersDirichlet;
  const auto invscheme   = pda::InviscidFluxReconstruction::Weno5;
  const auto visscheme   = pda::ViscousFluxReconstruction::FirstOrder;
  auto appObj      = pda::create_problem_eigen(meshObj, probid, invscheme, visscheme);
  using app_t = decltype(appObj);
  using state_t = typename app_t::state_type;
  using scalar_t = typename app_t::scalar_type;

  auto state = appObj.initialCondition();
  modify_state(state, meshObj);

  auto time = 0.0;
  auto velo = appObj.createVelocity();
  auto jac  = appObj.createJacobian();
  appObj.velocity(state, time, velo);
  appObj.jacobian(state, time, jac);

  writeToFileRank1(state, "state.txt");
  writeToFileRank1(velo,  "velo.txt");
  writeToFileSparseMat(jac,   "jacobian.txt");

  return 0;
}
