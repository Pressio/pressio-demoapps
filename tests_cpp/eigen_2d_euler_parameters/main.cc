
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
      std::cout << "*ERROR*: comparison failed for "
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
    const auto probId = pda::Euler2d::NormalShock;

    for (int api=0; api<4; ++api){
      if      (api==0){ std::cout << "\nNormalShock: testing overload: default params \n"; }
      else if (api==1){ std::cout << "\nNormalShock: testing overload: default params, customBCs \n"; }
      else if (api==2){ std::cout << "\nNormalShock: testing overload: default params, icFlag = 1 \n"; }
      else if (api==3){ std::cout << "\nNormalShock: testing overload: default params, icFlag = 1, customBCs \n"; }

      std::unordered_map<std::string, double> gold;
      gold["gamma"] = 1.4;
      gold["mach"] = 9.;

      if (api==0){
	auto appObj = pda::create_problem_eigen(meshObj, probId, flsc);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==1){
	auto noop = pda::impl::NoOperation<void>();
	auto appObj = pda::create_problem_eigen(meshObj, probId, flsc,
						noop, noop, noop, noop);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==2){
	auto appObj = pda::create_problem_eigen(meshObj, probId, flsc, 1);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==3){
	auto noop = pda::impl::NoOperation<void>();
	auto appObj = pda::create_problem_eigen(meshObj, probId, flsc,
						noop, noop, noop, noop, 1);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }
    }
  }

  {
    const auto probId = pda::Euler2d::NormalShock;
    const int icFlag = 1;
    for (int api=0; api<4; ++api){
      if      (api==0){ std::cout << "\nNormalShock: testing overload: custom params, icFlag =1, customBCs \n"; }

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
	  auto appObj = pda::create_problem_eigen(meshObj, probId, flsc,
						  noop, noop, noop, noop, icFlag, myparams);
	  if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
	}
      }
    }
  }

  {
    const auto probId = pda::Euler2d::Riemann;

    for (int api=0; api<6; ++api){
      if      (api==0){ std::cout << "\nRiemann: testing overload: default params \n"; }
      else if (api==1){ std::cout << "\nRiemann: testing overload: default params, customBCs \n"; }
      else if (api==2){ std::cout << "\nRiemann: testing overload: default params, icFlag = 1 \n"; }
      else if (api==3){ std::cout << "\nRiemann: testing overload: default params, icFlag = 1, customBCs \n"; }
      else if (api==4){ std::cout << "\nRiemann: testing overload: default params, icFlag = 2 \n"; }
      else if (api==5){ std::cout << "\nRiemann: testing overload: default params, icFlag = 2, customBCs \n"; }

      std::unordered_map<std::string, double> gold;
      gold["gamma"] = 1.4;

      if (api==0){
	gold["riemannTopRightPressure"] = 0.4;
	auto appObj = pda::create_problem_eigen(meshObj, probId, flsc);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==1){
	gold["riemannTopRightPressure"] = 0.4;
	auto noop = pda::impl::NoOperation<void>();
	auto appObj = pda::create_problem_eigen(meshObj, probId, flsc,
						noop, noop, noop, noop);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==2){
	gold["riemannTopRightPressure"] = 0.4;
	auto appObj = pda::create_problem_eigen(meshObj, probId, flsc, 1);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==3){
	gold["riemannTopRightPressure"] = 0.4;
	auto noop = pda::impl::NoOperation<void>();
	auto appObj = pda::create_problem_eigen(meshObj, probId, flsc,
						noop, noop, noop, noop, 1);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==4){
	gold["riemannTopRightPressure"] = 1.5;
	gold["riemannTopRightXVel"] = 0.0;
	gold["riemannTopRightYVel"] = 0.0;
	gold["riemannTopRightDensity"] = 1.5;
	gold["riemannBotLeftPressure"] = 0.029;
	auto appObj = pda::create_problem_eigen(meshObj, probId, flsc, 2);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }

      else if (api==5){
	gold["riemannTopRightPressure"] = 1.5;
	gold["riemannTopRightXVel"] = 0.0;
	gold["riemannTopRightYVel"] = 0.0;
	gold["riemannTopRightDensity"] = 1.5;
	gold["riemannBotLeftPressure"] = 0.029;
	auto noop = pda::impl::NoOperation<void>();
	auto appObj = pda::create_problem_eigen(meshObj, probId, flsc,
						noop, noop, noop, noop, 2);
	if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
      }
    }
  }

  {
    const auto probId = pda::Euler2d::Riemann;
    for (int icFlag = 1; icFlag <= 2; ++icFlag){
      for (int api=0; api<4; ++api){
	if (api==0){
	  std::cout << "\nRiemann: testing overload: custom params, icFlag = " << icFlag << " , customBCs \n";
	}

	// loop over each param and only change one at a time
	const int nparams = (icFlag == 1) ? 2 : 6;
	for (int i=0; i<nparams; ++i){
	  std::unordered_map<std::string, double> myparams;
	  if (i==0) myparams["gamma"] = 1.8;
	  if (i==1) myparams["riemannTopRightPressure"] = 8.3;
	  if (icFlag == 2) {
	    if (i==2) myparams["riemannTopRightXVel"] = 0.7;
	    if (i==3) myparams["riemannTopRightYVel"] = 1.2;
	    if (i==4) myparams["riemannTopRightDensity"] = 2.3;
	    if (i==5) myparams["riemannBotLeftPressure"] = 0.014;
	  }

	  std::unordered_map<std::string, double> gold;
	  gold["gamma"]  = (i==0) ? 1.8 : 1.4;
	  gold["riemannTopRightPressure"] = (i==1) ? 8.3 : (icFlag == 1) ? 0.4 : 1.5;
	  if (icFlag == 2) {
	    gold["riemannTopRightXVel"] = (i==2) ? 0.7 : 0.0;
	    gold["riemannTopRightYVel"] = (i==3) ? 1.2 : 0.0;
	    gold["riemannTopRightDensity"] = (i==4) ? 2.3 : 1.5;
	    gold["riemannBotLeftPressure"] = (i==5) ? 0.014 : 0.029;
	  }

	  if (api==0){
	    auto noop = pda::impl::NoOperation<void>();
	    auto appObj = pda::create_problem_eigen(meshObj, probId, flsc,
						    noop, noop, noop, noop, icFlag, myparams);
	    if (!check(appObj, gold)){ std::puts("FAILED"); return 0; }
	  }
	}
      }
    }
  }

  std::puts("PASS");
  return 0;
}
