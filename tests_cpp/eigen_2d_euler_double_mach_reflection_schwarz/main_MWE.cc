
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "../observer.hpp"

using namespace std;

int main()
{
    pressio::log::initialize(pressio::logto::terminal);
    pressio::log::setVerbosity({pressio::log::level::debug});

    namespace pda = pressiodemoapps;

    constexpr auto order   = pda::InviscidFluxReconstruction::FirstOrder;
    const auto probId  = pda::Euler2d::DoubleMachReflection;
    
    // ---- user inputs

    int ndomains = 2;
    int checkIdx = 0;

    // ---- end user inputs

    // getting types
    auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
    using mesh_t = decltype(meshObj);

    auto appObj = pda::create_problem_eigen(meshObj, probId, order);
    using app_t = decltype(appObj);
    using state_t = typename app_t::state_type;
    using jacob_t = typename app_t::jacobian_type;

    auto stepperObj = pressio::ode::create_implicit_stepper(
        pressio::ode::StepScheme::CrankNicolson, appObj);
    using lin_solver_t = pressio::linearsolvers::Solver<
    pressio::linearsolvers::iterative::Bicgstab, jacob_t>;
    lin_solver_t linSolverObj;
    auto NonLinSolver = pressio::nonlinearsolvers::create_newton_raphson(stepperObj, linSolverObj);
    NonLinSolver.setTolerance(1e-5);

    using nonlin_solver_t = decltype(NonLinSolver);
    using stepper_t = decltype(stepperObj);

    // problem objects for each domain
    vector<mesh_t> meshVec;
    vector<app_t> appVec;
    vector<state_t> stateVec;
    vector<stepper_t> stepperVec;
    vector<lin_solver_t> linSolverVec(ndomains);
    vector<nonlin_solver_t> nonlinSolverVec;
    for (int domIdx = 0; domIdx < ndomains; ++domIdx) {

        meshVec.push_back(pda::load_cellcentered_uniform_mesh_eigen("."));
        appVec.push_back(pda::create_problem_eigen(meshVec[domIdx], probId, order));
        stateVec.push_back(appVec[domIdx].initialCondition());
        stepperVec.push_back(pressio::ode::create_implicit_stepper(pressio::ode::StepScheme::CrankNicolson, appVec[domIdx]));
        nonlinSolverVec.push_back(pressio::nonlinearsolvers::create_newton_raphson(stepperVec[domIdx], linSolverVec[domIdx]));
        nonlinSolverVec[domIdx].setTolerance(1e-5);

    }

    const auto dtWrap = ::pressio::ode::StepSize<double>(0.0025);
    auto startTimeWrap = ::pressio::ode::StepStartAt<double>(0.01);
    const auto stepWrap = ::pressio::ode::StepCount(0);

    stepperVec[checkIdx](stateVec[checkIdx], startTimeWrap, stepWrap, dtWrap, nonlinSolverVec[checkIdx]);
    
    return 0;
}
