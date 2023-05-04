
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "../observer.hpp"

using namespace std;

int main()
{
  namespace pda  = pressiodemoapps;
  namespace plog = pressio::log;
  namespace pode = pressio::ode;
  namespace pls  = pressio::linearsolvers;
  namespace pnls = pressio::nonlinearsolvers;

  plog::initialize(pressio::logto::terminal);
  plog::setVerbosity({pressio::log::level::debug});

  // ---- user inputs
  int ndomains = 2;
  const auto order  = pda::InviscidFluxReconstruction::FirstOrder;
  const auto probId = pda::Euler2d::DoubleMachReflection;
  const auto scheme = pode::StepScheme::CrankNicolson;
  // ---- end user inputs

  // aliases
  using mesh_t  = decltype( pda::load_cellcentered_uniform_mesh_eigen(".") );
  using app_t   = decltype( pda::create_problem_eigen(std::declval<mesh_t>(), probId, order) );
  using state_t = typename app_t::state_type;
  using jacob_t = typename app_t::jacobian_type;
  using lin_solver_t = pls::Solver<pls::iterative::Bicgstab, jacob_t>;
  using stepper_t = decltype( pode::create_implicit_stepper(scheme, std::declval<app_t &>()) );
  using nonlin_solver_t = decltype( pnls::create_newton_raphson( std::declval<stepper_t &>(), std::declval<lin_solver_t&>()) );

  //
  // create things
  //
  /* Note:
     only need a single linear solver instance, all the nonlinear solvers will reference the same
     but make sure it does not go out of scope */
  lin_solver_t linSolverObj;
  std::string path = ".";
  std::vector<mesh_t>     meshVec(ndomains);
  std::vector<app_t>      appVec(ndomains);
  std::vector<state_t>    stateVec(ndomains);
  std::vector<stepper_t>  stepperVec;
  vector<nonlin_solver_t> nonlinSolverVec;
  std::cout << "INIT \n";
  for (int domIdx = 0; domIdx < ndomains; ++domIdx)
  {
    meshVec[domIdx]  = pda::load_cellcentered_uniform_mesh_eigen(path);
    appVec[domIdx]   = pda::create_problem_eigen(meshVec[domIdx], probId, order);
    stateVec[domIdx] = appVec[domIdx].initialCondition();

    stepperVec.emplace_back( pode::create_implicit_stepper(scheme, appVec[domIdx]) );
    nonlinSolverVec.emplace_back( pnls::create_newton_raphson(stepperVec.back(), linSolverObj) );

    std::cout << " LOOP: "
	      << "domIdx = " << domIdx << " "
	      << "&meshObj = " << &meshVec[domIdx] << " "
	      << "&appObj =  " << &appVec[domIdx]
	      << std::endl;
  }

  std::cout << "\nVERIFY \n";
  for (int domIdx = 0; domIdx < ndomains; ++domIdx) {  //
    std::cout << " LOOP: "
	      << "domIdx = " << domIdx << " "
	      << "&meshObj = " << &meshVec[domIdx] << " "
	      << "&appObj =  " << &appVec[domIdx]
	      << std::endl;
  }

  //
  // try solve
  //
  std::cout << "\n";
  auto fixTol   = [](nonlin_solver_t & s){ s.setTolerance(1e-5); };
  auto fixIters = [](nonlin_solver_t & s){ s.setMaxIterations(3); };
  std::for_each(nonlinSolverVec.begin(), nonlinSolverVec.end(), fixIters);

  int numSteps = 3;
  const auto dtWrap = pode::StepSize<double>(0.0025);
  const auto startTimeWrap = pode::StepStartAt<double>(0);
  double time = startTimeWrap.get();
  for (int stepId=1; stepId<=numSteps; ++stepId)
    {
      std::cout << " Doing step = " << stepId << "\n";
      for (int domIdx = 0; domIdx < ndomains; ++domIdx) {
	std::cout << " Doing domain = " << domIdx << "\n";
	stepperVec[domIdx](stateVec[domIdx], startTimeWrap,
			   pode::StepCount(stepId), dtWrap,
			   nonlinSolverVec[domIdx]);
      }
      time += dtWrap.get();
    }

  return 0;
}
