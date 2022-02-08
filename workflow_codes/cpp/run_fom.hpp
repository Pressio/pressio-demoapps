
#ifndef RUN_FOM_HPP_
#define RUN_FOM_HPP_

#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "observer.hpp"
#include "rom_shared.hpp"

template<class FomType, class ParserType>
void run_fom(FomType & appObj, const ParserType & parser)
{
  namespace pode = pressio::ode;

  // initial condition
  using state_t = typename FomType::state_type;
  state_t state = appObj.initialCondition();

  StateObserver<state_t> stateObs("fom_snaps_state.bin", parser.fomStateSamplingFreq());
  const auto startTime = static_cast<typename FomType::scalar_type>(0);

  const auto odeScheme = parser.fomOdeScheme();
  if (pressio::ode::is_explicit_scheme(odeScheme))
  {
    // observer to sample rhs
    VelocityObserver<> veloObs("fom_snaps_rhs.bin", parser.fomVelocitySamplingFreq());

    auto stepperObj = pode::create_explicit_stepper(odeScheme, state, appObj);
    pode::advance_n_steps_and_observe(stepperObj, state, startTime,
				      parser.fomTimeStepSize(),
				      parser.fomNumSteps(),
				      stateObs, veloObs);
  }
  else{
    namespace plins  = pressio::linearsolvers;
    namespace pnlins = pressio::nonlinearsolvers;

    auto stepperObj = pode::create_implicit_stepper(odeScheme, state, appObj);

    using jacob_t      = typename FomType::jacobian_type;
    using lin_solver_t = plins::Solver<plins::iterative::Bicgstab, jacob_t>;
    lin_solver_t linSolverObj;

    // make nonlinear solver
    const auto nonlinSolverType = parser.fomNonlinearSolverType();
    if (nonlinSolverType != "NewtonRaphson"){
      throw std::runtime_error("FOM requires NewtonRaphson");
    }

    auto NonLinSolver = pnlins::create_newton_raphson(stepperObj, state, linSolverObj);
    const auto tol = parser.fomNonlinearSolverTolerance();
    NonLinSolver.setTolerance(tol);

    pode::advance_n_steps_and_observe(stepperObj, state, startTime,
				      parser.fomTimeStepSize(),
				      parser.fomNumSteps(),
				      stateObs, NonLinSolver);
  }
}

#endif
