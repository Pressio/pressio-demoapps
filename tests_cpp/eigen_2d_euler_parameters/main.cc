
#include "pressiodemoapps/euler2d.hpp"

namespace
{
template<class AppT>
bool check(const AppT & app,
	   const std::unordered_map<std::string, double> & gold)
{

  for (auto it=gold.cbegin(); it!=gold.cend(); ++it){
    const auto value = app.queryParameter(it->first);
    const auto diff = std::abs(value - it->second);
    if ( diff > 1e-13 ){
      std::cout << "*ERROR*: comparison failed for"
		<< it->first << " : "
		<< " found = " << value << " gold  = " << it->second
		<< '\n';
      return false;
    }
  }
  return true;
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
    for (int api=0; api<4; ++api){
      if      (api==0){ std::cout << "\ntesting overload: default params \n"; }
      else if (api==1){ std::cout << "\ntesting overload: default params, customBCs \n"; }
      else if (api==2){ std::cout << "\ntesting overload: default params, icFlag = 1 \n"; }
      else if (api==3){ std::cout << "\ntesting overload: default params, icFlag = 1, customBCs \n"; }

      std::unordered_map<std::string, double> gold;
      gold["gamma"] = 1.4;
      gold["mach"] = 9.;

      if (api==0){
	auto appObj = pda::create_problem_eigen(meshObj, pda::Euler2d::NormalShock, flsc);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==1){
	auto noop = pda::impl::NoOperation<void>();
	auto appObj = pda::create_problem_eigen(meshObj, pda::Euler2d::NormalShock, flsc,
						noop, noop, noop, noop);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==2){
	auto appObj = pda::create_problem_eigen(meshObj, pda::Euler2d::NormalShock, flsc, 1);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==3){
	auto noop = pda::impl::NoOperation<void>();
	auto appObj = pda::create_problem_eigen(meshObj, pda::Euler2d::NormalShock, flsc,
						noop, noop, noop, noop, 1);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }
    }
  }

  {
    const int icFlag = 1;
    for (int api=0; api<4; ++api){
      if      (api==0){ std::cout << "\ntesting overload: custom params, icFlag =1, customBCs \n"; }

      // loop over each param and only change one at a time
      for (int i=0; i<2; ++i)
      {
	std::unordered_map<std::string, double> myparams;
	if (i==0) myparams["gamma"] = 1.8;
	if (i==1) myparams["mach"]  = 12.;

	std::unordered_map<std::string, double> gold;
	gold["gamma"]  = (i==0) ? 1.8 : 1.4;
	gold["mach"] = (i==1) ? 12. : 9.0;

	if (api==0){
	  auto noop = pda::impl::NoOperation<void>();
	  auto appObj = pda::create_problem_eigen(meshObj, pda::Euler2d::NormalShock, flsc,
						  noop, noop, noop, noop, icFlag, myparams);
	  if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
	}
      }
    }
  }

  std::puts("PASS");
  return 0;
}
