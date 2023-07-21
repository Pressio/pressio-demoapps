#include <vector>
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/advection1d.hpp"
#include "../observer.hpp"

template <class AppT, class MeshT, class ...AppArgs>
auto create_apps(std::vector<MeshT> const & meshes, AppArgs && ... appargs)
{
  namespace pda  = pressiodemoapps;
  std::vector<AppT> r;
  for (std::size_t i=0; i<meshes.size(); ++i){
    r.emplace_back( pda::create_problem_eigen(meshes[i], std::forward<AppArgs>(appargs)...) );
  }
  return r;
}

template <class MeshT>
auto create_meshes(int n, std::string const & meshRoot)
{
  namespace pda = pressiodemoapps;
  std::vector<MeshT> r;
  for (int i=0; i<n; ++i){
    r.emplace_back( pda::load_cellcentered_uniform_mesh_eigen(meshRoot + "/domain" + std::to_string(i)) );
  }
  return r;
}

template<class AppType, class MeshT>
class SchwarzDecomp
{
  using mesh_t    = MeshT;
  using state_t   = typename AppType::state_type;
  using jacob_t   = typename AppType::jacobian_type;
  using stepper_t = decltype( pressio::ode::create_implicit_stepper(pressio::ode::StepScheme(),
								    std::declval<AppType &>()) );
  using lin_solver_tag = pressio::linearsolvers::iterative::Bicgstab;
  using lin_solver_t   = pressio::linearsolvers::Solver<lin_solver_tag, jacob_t>;
  using nonlin_solver_t = decltype( pressio::create_newton_solver( std::declval<stepper_t &>(),
								   std::declval<lin_solver_t&>()) );
public:
  template<class ...Args>
  SchwarzDecomp(int n,
		const std::string & meshRoot,
		pressio::ode::StepScheme scheme,
		Args && ...args)
    : ndomains(n)
    , meshVec( create_meshes<mesh_t>(n, ".") )
    , appVec( create_apps<AppType>(meshVec, std::forward<Args>(args)...) )
  {
    for (int domIdx = 0; domIdx < n; ++domIdx){
      stateVec.emplace_back( appVec[domIdx].initialCondition() );
      stepperVec.emplace_back( pressio::ode::create_implicit_stepper(scheme, appVec[domIdx]) );
      nonlinSolverVec.emplace_back( pressio::create_newton_solver(stepperVec.back(), linSolverObj) );
    }
  }

  void dummy(){
    namespace pode = pressio::ode;
    auto fixIters = [](nonlin_solver_t & s){ s.setMaxIterations(3); };
    std::for_each(nonlinSolverVec.begin(), nonlinSolverVec.end(), fixIters);

    int numSteps = 3;
    const auto dtWrap = pode::StepSize<double>(0.0025);
    const auto startTimeWrap = pode::StepStartAt<double>(0);
    for (int stepId=1; stepId<=numSteps; ++stepId)
      {
	std::cout << " Doing step = " << stepId << "\n";
	for (int domIdx = 0; domIdx < ndomains; ++domIdx) {
	  std::cout << " Doing domain = " << domIdx << "\n";
	  stepperVec[domIdx](stateVec[domIdx], startTimeWrap,
			     pode::StepCount(stepId), dtWrap,
			     nonlinSolverVec[domIdx]);
	}
      }
  }

private:
  int ndomains;
  lin_solver_t linSolverObj = {};
  std::vector<mesh_t> meshVec;
  std::vector<AppType> appVec;
  std::vector<state_t> stateVec;
  std::vector<stepper_t> stepperVec;
  std::vector<nonlin_solver_t> nonlinSolverVec;
  //std::vector<std::vector<state_t>> stateHistVec;
};

int main()
{
  namespace pda  = pressiodemoapps;
  namespace plog = pressio::log;
  namespace pode = pressio::ode;
  namespace pls  = pressio::linearsolvers;
  namespace pnls = pressio::nonlinearsolvers;

  plog::initialize(pressio::logto::terminal);
  plog::setVerbosity({pressio::log::level::debug});

  int ndomains = 3;
  const auto order  = pda::InviscidFluxReconstruction::FirstOrder;
  const auto probId = pda::Advection1d::PeriodicLinear;
  const auto scheme = pode::StepScheme::CrankNicolson;
  using mesh_t = pressiodemoapps::cellcentered_uniform_mesh_eigen_type;
  using app_t = decltype( pda::create_problem_eigen(std::declval<mesh_t>(), probId, order) );

  using decomp_type = SchwarzDecomp<app_t, mesh_t>;
  decomp_type d(ndomains, ".", scheme, probId, order);

  // this test passes if things "run" and nothing is thrown
  try{
    d.dummy();
  }
  catch (...){
    std::puts("FAILED");
  }
  std::puts("PASS");

  return 0;
}
