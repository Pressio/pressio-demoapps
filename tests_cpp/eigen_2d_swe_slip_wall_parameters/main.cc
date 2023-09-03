
#include "pressiodemoapps/swe2d.hpp"

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
    for (int api=0; api<6; ++api){
      if      (api==0){ std::cout << "\ntesting overload: default params \n"; }
      else if (api==1){ std::cout << "\ntesting overload: default params, customBCs \n"; }
      else if (api==2){ std::cout << "\ntesting overload: default params, icFlag = 1 \n"; }
      else if (api==3){ std::cout << "\ntesting overload: default params, icFlag = 1, customBCs \n"; }
      else if (api==4){ std::cout << "\ntesting overload: default params, icFlag = 2 \n"; }
      else if (api==5){ std::cout << "\ntesting overload: default params, icFlag = 2, customBCs \n"; }

      std::unordered_map<std::string, double> gold;
      if (api<4){
	gold["gravity"]	       = 9.8;
	gold["coriolis"]       = -3;
	gold["pulseMagnitude"] = 1./8.;
	gold["pulseX"]	       = 1.0;
	gold["pulseY"]	       = 1.0;
      }
      else{
	gold["gravity"]	        = 9.8;
	gold["coriolis"]        = -3.;
	gold["pulseMagnitude1"] = 1./10.;
	gold["pulseX1"]	        = -2.;
	gold["pulseY1"]	        = -2.;
	gold["pulseMagnitude2"] = 1./8.;
	gold["pulseX2"]	        = 2.;
	gold["pulseY2"]         = 2.;
      }

      if (api==0){
	auto appObj = pda::create_problem_eigen(meshObj, pda::Swe2d::SlipWall, flsc);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==1){
	auto noop = pda::impl::NoOperation<void>();
	auto appObj = pda::create_problem_eigen(meshObj, pda::Swe2d::CustomBCs, flsc,
						noop, noop, noop, noop);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==2){
	auto appObj = pda::create_problem_eigen(meshObj, pda::Swe2d::SlipWall, flsc, 1);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==3){
	auto noop = pda::impl::NoOperation<void>();
	auto appObj = pda::create_problem_eigen(meshObj, pda::Swe2d::CustomBCs, flsc,
						noop, noop, noop, noop, 1);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==4){
	auto appObj = pda::create_problem_eigen(meshObj, pda::Swe2d::SlipWall, flsc, 2);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==5){
	auto noop = pda::impl::NoOperation<void>();
	auto appObj = pda::create_problem_eigen(meshObj, pda::Swe2d::CustomBCs, flsc,
						noop, noop, noop, noop, 2);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }
    }
  }

  {
    const int icFlag = 1;
    for (int api=0; api<4; ++api){
      if      (api==0){ std::cout << "\ntesting overload: custom params, icFlag =1 \n"; }
      else if (api==1){ std::cout << "\ntesting overload: custom params, icFlag =1, customBCs \n"; }

      // loop over each param and only change one at a time
      for (int i=0; i<5; ++i)
      {
	std::unordered_map<std::string, double> myparams;
	if (i==0) myparams["gravity"]        = 10.;
	if (i==1) myparams["coriolis"]	   = -5.;
	if (i==2) myparams["pulseMagnitude"] = 0.05;
	if (i==3) myparams["pulseX"]	   = -1.5;
	if (i==4) myparams["pulseY"]	   = -1.8;

	std::unordered_map<std::string, double> gold;
	gold["gravity"]        = (i==0) ? 10   : 9.8;
	gold["coriolis"]       = (i==1) ? -5.  : -3.0;
	gold["pulseMagnitude"] = (i==2) ? 0.05 : 1./8.;
	gold["pulseX"]	     = (i==3) ? -1.5 : 1.;
	gold["pulseY"]	     = (i==4) ? -1.8 : 1.;

	if (api==0){
	  auto appObj = pda::create_problem_eigen(meshObj, pda::Swe2d::SlipWall, flsc, icFlag, myparams);
	  if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
	}
	else if (api==1){
	  auto noop = pda::impl::NoOperation<void>();
	  auto appObj = pda::create_problem_eigen(meshObj, pda::Swe2d::CustomBCs, flsc,
						  noop, noop, noop, noop, icFlag, myparams);
	  if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
	}
      }
    }
  }

  {
    const int icFlag = 2;
    for (int api=0; api<4; ++api){
      if      (api==0){ std::cout << "\ntesting overload: custom params, icFlag =2\n"; }
      else if (api==1){ std::cout << "\ntesting overload: custom params, icFlag =2, customBCs \n"; }

      // loop over each param and only change one at a time
      for (int i=0; i<8; ++i)
      {
	std::unordered_map<std::string, double> myparams;
	if (i==0) myparams["gravity"]         = 10.;
	if (i==1) myparams["coriolis"]	      = -5.;
	if (i==2) myparams["pulseMagnitude1"] = 0.05;
	if (i==3) myparams["pulseX1"]	      = -1.5;
	if (i==4) myparams["pulseY1"]	      = -1.8;
	if (i==5) myparams["pulseMagnitude2"] = 0.065;
	if (i==6) myparams["pulseX2"]	      = -2.5;
	if (i==7) myparams["pulseY2"]	      = -2.8;

	std::unordered_map<std::string, double> gold;
	gold["gravity"]         = (i==0) ? 10   : 9.8;
	gold["coriolis"]        = (i==1) ? -5.  : -3.0;
	gold["pulseMagnitude1"] = (i==2) ? 0.05 : 1./10.;
	gold["pulseX1"]	        = (i==3) ? -1.5 : -2.;
	gold["pulseY1"]	        = (i==4) ? -1.8 : -2.;
	gold["pulseMagnitude2"] = (i==5) ? 0.065 : 1./8.;
	gold["pulseX2"]	        = (i==6) ? -2.5 : 2.;
	gold["pulseY2"]	        = (i==7) ? -2.8 : 2.;

	if (api==0){
	  auto appObj = pda::create_problem_eigen(meshObj, pda::Swe2d::SlipWall, flsc, icFlag, myparams);
	  if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
	}
	else if (api==1){
	  auto noop = pda::impl::NoOperation<void>();
	  auto appObj = pda::create_problem_eigen(meshObj, pda::Swe2d::CustomBCs, flsc,
						  noop, noop, noop, noop, icFlag, myparams);
	  if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
	}
      }
    }
  }

  std::puts("PASS");
  return 0;
}
