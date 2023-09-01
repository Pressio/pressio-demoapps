
#include "pressiodemoapps/swe2d.hpp"

namespace
{
template<class AppT>
void check(const AppT & app,
	   const std::unordered_map<std::string, double> & gold)
{

  for (auto it=gold.cbegin(); it!=gold.cend(); ++it){
    const auto value = app.queryParameter(it->first);
    const auto diff = std::abs(value - it->second);
    std::cout << it->first
	      << " : "
	      << " found = " << value
	      << " gold  = " << it->second
	      << '\n';
    if ( diff > 1e-13 ){
      const std::string msg = "comparison for: " + it->first + ", failed";
      throw std::runtime_error(msg);
    }
  }
}

}

int main()
{
  namespace pda = pressiodemoapps;
  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
  // this test only cares about the physical parameters,
  // it does not matter here which flux scheme we use
  const auto flsc = pda::InviscidFluxReconstruction::FirstOrder;

  {
    std::cout << "\ntesting overload 1\n";
    const auto probId  = pda::Swe2d::SlipWall;
    auto appObj = pda::create_problem_eigen(meshObj, probId, flsc);
    std::unordered_map<std::string, double> gold;
    gold["gravity"] = 9.8;
    gold["coriolis"] = static_cast<double>(-3);
    gold["pulseMagnitude"] = static_cast<double>(1)/8;
    gold["pulseX"] = static_cast<double>(1);
    gold["pulseY"] = static_cast<double>(1);
    check(appObj, gold);
  }

  {
    std::cout << "\ntesting overload 2\n";
    auto appObj = pda::  create_slip_wall_swe_2d_problem_eigen(meshObj, flsc,
							       12., 4.4, 0.5);
    std::unordered_map<std::string, double> gold;
    gold["gravity"] = 12.;
    gold["coriolis"] = static_cast<double>(4.4);
    gold["pulseMagnitude"] = static_cast<double>(1)/2;
    gold["pulseX"] = static_cast<double>(1);
    gold["pulseY"] = static_cast<double>(1);
    check(appObj, gold);
  }

  {
    std::cout << "\ntesting overload 3\n";
    const auto probId  = pda::Swe2d::CustomBCs;
    auto noop = pda::impl::NoOperation<void>();
    auto appObj = pda::create_problem_eigen(meshObj, probId, flsc,
					    noop, noop, noop, noop);
    std::unordered_map<std::string, double> gold;
    gold["gravity"] = 9.8;
    gold["coriolis"] = static_cast<double>(-3);
    gold["pulseMagnitude"] = static_cast<double>(1)/8;
    gold["pulseX"] = static_cast<double>(1);
    gold["pulseY"] = static_cast<double>(1);
    check(appObj, gold);
  }

  {
    std::cout << "\ntesting overload 4\n";
    const auto probId  = pda::Swe2d::CustomBCs;
    auto noop = pda::impl::NoOperation<void>();
    std::unordered_map<std::string, double> myparams;
    myparams["gravity"] = 10;
    myparams["coriolis"] = static_cast<double>(-5);
    myparams["pulseMagnitude"] = static_cast<double>(1.5)/8;
    myparams["pulseX"] = static_cast<double>(1.5);
    myparams["pulseY"] = static_cast<double>(1.5);

    auto appObj = pda::create_problem_eigen(meshObj, probId, flsc,
					    noop, noop, noop, noop,
					    myparams);
    check(appObj, myparams);
  }


  return 0;
}
