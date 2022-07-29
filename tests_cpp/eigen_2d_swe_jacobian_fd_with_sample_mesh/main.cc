
#include "pressio/type_traits.hpp"
#include "pressiodemoapps/swe2d.hpp"
#include <iomanip>
#include <random>

template<class T>
void writeToFile(const T& obj, const std::string & fileName)
{
  std::ofstream file; file.open(fileName);
  for (int i=0; i<obj.size(); i++){
    file << std::setprecision(16) << obj(i) << " \n";
  }
  file.close();
}

struct UnifDist{
  using dist_type = std::uniform_real_distribution<double>;
  std::random_device rd;
  std::mt19937 m_gen{rd()};
  dist_type m_dist;

  UnifDist(const double a, const double b) : m_dist(a, b){}
  double operator()() { return m_dist(m_gen); }
};

template<class state_type, class mesh_t>
void modify_state(state_type & state,
      const mesh_t & meshObj)
{
  UnifDist randObj(0.1, 0.9);

  const auto & x = meshObj.viewX();
  for (int i=0; i<::pressiodemoapps::extent(x,0); ++i){
    const auto ind = i*3;
    // const auto pert = 20.*std::sin(x(i)*y(i));
    // std::cout << x(i) << " " << y(i) << " " << pert << '\n';
    state(ind)   = 1.0 + randObj();
    state(ind+1) = 1.5 + randObj();
    state(ind+2) = 2.0 + randObj();
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
#else
  const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probId  = pda::Swe2d::SlipWall;
  auto appObj = pda::create_problem_eigen(meshObj, probId, order);

  auto state = appObj.initialCondition();
  modify_state(state, meshObj);
  writeToFile(state, "IC.txt");

  const double eps = 1e-8;
  auto velo = appObj.createRightHandSide();
  auto J = appObj.createJacobian();

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
    decltype(velo) velo2(velo.size());
    appObj.rightHandSide(state2, 0., velo2);

    auto Ja_fd = (velo2 - velo)/eps;
    for (int i=0; i<Ja.size(); ++i)
    {
      const auto diff = std::abs(Ja(i)- Ja_fd(i));
      if (loop == 0){
	       printf(" i=%2d J(i)=%16.12f J_fd(i)=%16.12f diff=%e \n",
	       i, Ja(i), Ja_fd(i), diff);
      }

      if (diff > 1e-5){
	       std::puts("FAILED");
	       //return 0;
      }
    }

    // second order
    auto state3 = state-eps*a;
    decltype(velo) velo3(velo.size());
    appObj.rightHandSide(state3, 0., velo3);
    auto Ja_fd_2 = (velo2 - velo3)/(2.*eps);
    for (int i=0; i<Ja.size(); ++i){
      const auto diff = std::abs(Ja(i)- Ja_fd_2(i));
      if (loop == 0){
    	printf(" i=%2d J(i)=%10.6f J_fd(i)=%10.6f diff=%e \n", i, Ja(i), Ja_fd(i), diff);
      }

      if (diff > 1e-5){
    	  std::puts("FAILED");
    	  return 0;
      }
    }

  }

  std::puts("PASS");
  return 0;
}
