
#include "pressiodemoapps/euler3d.hpp"
#include <iomanip>
#include <array>

template<class T>
void writeToFile(const T& obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (int i=0; i<obj.size(); i++){
    file << std::setprecision(14) << obj(i) << " \n";
  }
  file.close();
}

// template<class state_type, class mesh_t>
// void modify_state(state_type & state,
// 		  const mesh_t & meshObj)
// {

//   using scalar_type = typename mesh_t::scalar_t;
//   const auto gamma = 1.4;
//   constexpr int numDofPerCell = 5;
//   constexpr auto zero = static_cast<scalar_type>(0);
//   constexpr auto one  = static_cast<scalar_type>(1);
//   constexpr auto two = static_cast<scalar_type>(2);

//   const auto & x= meshObj.viewX();
//   const auto & y= meshObj.viewY();
//   const auto & z= meshObj.viewZ();
//   const auto gammaMinusOne = gamma - one;
//   const auto gammaMinusOneInv = one/gammaMinusOne;

//   std::array<scalar_type, numDofPerCell> prim;
//   for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
//     {
//       const auto sx = std::sin(8.*M_PI*x(i));
//       const auto sy = std::sin(12.*M_PI*y(i));
//       const auto sz = std::sin(6.*M_PI*z(i));
//       const auto pert = 0.001*(sx+sy+sz);
//       prim[0] += pert;
//       prim[1] += pert;
//       prim[2] += pert;
//       prim[3] += pert;

//       const auto ind = i*numDofPerCell;
//       state(ind)   = prim[0];
//       state(ind+1) = prim[0]*prim[1];
//       state(ind+2) = prim[0]*prim[2];
//       state(ind+3) = pressiodemoapps::eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv, prim);
//     }
// }

int main()
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");

#if defined USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#elif defined USE_FIRSTORDER
  constexpr auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  auto appObj = pda::create_problem_eigen(meshObj, pda::Euler3d::PeriodicSmooth, order);
  using app_t = decltype(appObj);
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;

  app_state_t state = appObj.initialCondition();
  //modify_state(state, meshObj);
  writeToFile(state, "IC.txt");

  auto velo = appObj.createRightHandSide();
  auto J = appObj.createJacobian();
  const double eps = 1e-8;

  // make sure repeated evaluations work
  // not just a single time
  for (int loop=0; loop<5; ++loop)
  {

    appObj.rightHandSide(state, 0., velo);
    appObj.jacobian(state, 0., J);

     Eigen::VectorXd a = Eigen::VectorXd::Random(state.size());
     auto Ja = J*a;

    // first order
    auto state2 = state+eps*a;
    app_rhs_t velo2(velo.size());
    appObj.rightHandSide(state2, 0., velo2);

    app_rhs_t Ja_fd = (velo2 - velo)/eps;
    for (int i=0; i<Ja.size(); ++i)
      {
	const auto diff = std::abs(Ja(i)- Ja_fd(i));
	if (loop==1){
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
    appObj.rightHandSide(state3, 0., velo3);
    auto Ja_fd_2 = (velo2 - velo3)/(2.*eps);
    for (int i=0; i<Ja.size(); ++i){
      const auto diff = std::abs(Ja(i)- Ja_fd_2(i));
      if (diff > 1e-5){
        std::puts("FAILED");
        return 0;
      }
    }

  }

  std::puts("PASS");
  return 0;
}
