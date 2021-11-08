
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

int main(int argc, char *argv[])
{

  namespace pda = pressiodemoapps;
  const auto meshObj = pda::loadCellCenterUniformMeshEigen(".");
  constexpr auto order   = pda::InviscidFluxReconstruction::FirstOrder;

  // here we want to test FD Jacobian when we have Neumann BC on x
  // and reflective on y
  // we can use NormalShock for this since that problem has the correct BCs

  auto appObj = pda::createProblemEigen(meshObj, pda::Euler2d::NormalShock, order);
  using app_t = decltype(appObj);
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  using app_jacob_t	= typename app_t::jacobian_type;

  app_state_t state = appObj.initialCondition();
  state.setZero();
  const auto & x= meshObj.viewX();
  const auto & y= meshObj.viewY();
  std::array<scalar_t, 4> prim;
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
    {
      const auto pert = 0.01*std::sin(16.*M_PI*x(i))*std::sin(8.*M_PI*y(i));

      if (x(i) <= 1.0){
	prim[0] = 1;
	prim[1] = 0.1;
	prim[2] = 0.;
	prim[3] = 1;
      }

      else{
	prim[0] = static_cast<scalar_t>(0.125);
	prim[1] = 0.1;
	prim[2] = 0.;
	prim[3] = static_cast<scalar_t>(0.1);
      }

      prim[0] += pert;
      prim[1] += pert;
      prim[2] += pert;
      prim[3] += pert;

      const auto ind = i*4;
      state(ind)   = prim[0];
      state(ind+1) = prim[0]*prim[1];
      state(ind+2) = prim[0]*prim[2];
      state(ind+3) = pressiodemoapps::eulerEquationsComputeEnergyFromPrimitive(1.4, prim);
    }
  writeToFile(state, "IC.txt");

  const scalar_t eps = 1e-8;

  auto velo = appObj.createVelocity();
  auto J = appObj.createJacobian();
  appObj.velocity(state, 0., velo);
  appObj.jacobian(state, 0., J);
  Eigen::VectorXd a = Eigen::VectorXd::Random(state.size());

  // first order
  Eigen::VectorXd state2(state.size());
  for (int i=0; i<state2.size(); ++i){
    state2(i) = state(i) + eps * a(i);
  }

  app_rhs_t velo2(velo.size());
  appObj.velocity(state2, 0., velo2);

  auto Ja = J*a;
  Eigen::VectorXd Ja_fd(Ja.size());
  for (int i=0; i<Ja_fd.size(); ++i){
    Ja_fd(i) = (velo2(i) - velo(i))/eps;
  }
  for (int i=0; i<Ja.size(); ++i)
  {
    const auto diff = std::abs(Ja(i)- Ja_fd(i));
    std::cout << i << " "
	      << Ja(i) << " "
	      << Ja_fd(i) << " "
	      << diff << " "
	      << std::endl;

    if ((i+1) % 4 ==0){
      std::cout << "\n";
    }

    if (diff > 1e-4){
      std::puts("FAILED");
      return 0;
    }
  }

  // // second order
  // auto state3 = state-eps*a;
  // app_rhs_t velo3(velo.size());
  // appObj.velocity(state3, 0., velo3);
  // auto Ja_fd_2 = (velo2 - velo3)/(2.*eps);
  // for (int i=0; i<Ja.size(); ++i){
  //   const auto diff = std::abs(Ja(i)- Ja_fd_2(i));
  //   if (diff > 1e-6){
  //     std::puts("FAILED");
  //     return 0;
  //   }
  // }

  std::puts("PASS");
  return 0;
}
