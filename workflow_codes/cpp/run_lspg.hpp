
#ifndef RUN_LSPG_HPP_
#define RUN_LSPG_HPP_

#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressio/rom_lspg.hpp"
#include "observer.hpp"
#include "rom_shared.hpp"

template<class IntType>
IntType find_index(const Eigen::Matrix<IntType, -1, 1> & vector,
		   IntType target)
{
  for (IntType j=0; j<vector.size(); ++j){
    if (vector[j] == target){
      return j;
    }
  }

  return std::numeric_limits<IntType>::max();
}

template<class ScalarType>
struct HypRedUpdater
{
  using vec_operand_type = Eigen::Matrix<ScalarType, -1, 1>;
  using mat_ll_operand_type = Eigen::Matrix<ScalarType, -1, -1>;
  std::vector<int> indices_;
  int numDofsPerCell_ = {};

  template<class SampleMeshGidsType, class StencilMeshGidsType>
  HypRedUpdater(const SampleMeshGidsType & smGids,
		const StencilMeshGidsType & stencilGids,
		int numDofsPerCell)
    : indices_(smGids.size()), numDofsPerCell_(numDofsPerCell)
  {
    std::cout << "updater \n";
    for (int i=0; i<indices_.size(); ++i){
      const auto index = find_index<int>(stencilGids, smGids[i]);
      assert(index != std::numeric_limits<int>::max());
      indices_[i] = index;
    }

    // std::ofstream file; file.open("indices.txt");
    // for (std::size_t i=0; i<indices_.size(); i++){
    //   file << std::setprecision(14) << indices_[i] << " \n";
    // }
    // file.close();
  }

  // a = alpha*a + beta*b (a,b potentially with different distributions)
  void updateSampleMeshOperandWithStencilMeshOne(vec_operand_type & a,
						 ScalarType alpha,
						 const vec_operand_type & b,
						 ScalarType beta) const
  {
    for (std::size_t i=0; i<indices_.size(); ++i){
      const std::size_t r = i*numDofsPerCell_;
      const std::size_t g = indices_[i]*numDofsPerCell_;
      for (std::size_t k=0; k<numDofsPerCell_; ++k){
	a(r+k) = alpha*a(r+k) + beta*b(g+k);
      }
    }
  }

  void updateSampleMeshOperandWithStencilMeshOne(mat_ll_operand_type & a,
						 ScalarType alpha,
						 const mat_ll_operand_type & b,
						 ScalarType beta) const
  {
    for (std::size_t j=0; j<b.cols(); ++j)
    {
      for (std::size_t i=0; i<indices_.size(); ++i)
      {
	const std::size_t r = i*numDofsPerCell_;
	const std::size_t g = indices_[i]*numDofsPerCell_;
	for (std::size_t k=0; k<numDofsPerCell_; ++k)
	{
	  a(r+k,j) = alpha*a(r+k,j) + beta*b(g+k,j);
	}
      }
    }
  }
};

template<
  class RomProblemType,
  class ParserType,
  class RomStateType
  >
void pick_solver_and_run(const ParserType & parser,
			 RomProblemType & problem,
			 RomStateType & romState)

{
  using scalar_t = typename RomProblemType::scalar_type;
  namespace pode = pressio::ode;
  namespace pnlins = pressio::nonlinearsolvers;

  // create observer
  StateObserver<RomStateType> obs("rom_snaps.bin", parser.romStateSamplingFreq());
  const int nSteps = parser.romNumSteps();

  if (parser.romNonlinearSolverType() == "GaussNewton")
  {
    using hessian_t = Eigen::Matrix<scalar_t, -1, -1>;
    // using solver_tag = pressio::linearsolvers::iterative::Bicgstab;
    using solver_tag = pressio::linearsolvers::direct::HouseholderQR;
    using linear_solver_t = pressio::linearsolvers::Solver<solver_tag, hessian_t>;
    linear_solver_t linSolverObj;

    auto solver = pnlins::create_gauss_newton(problem, romState, linSolverObj);
    solver.setStoppingCriterion(pnlins::Stop::WhenGradientAbsoluteNormBelowTolerance);

    const auto tol = parser.romNonlinearSolverTolerance();
    solver.setTolerance(tol);
    //solver.setUpdatingCriterion(pnlins::Update::Armijo);

    pode::advance_n_steps_and_observe(problem, romState, 0.,
				      parser.romTimeStepSize(),
				      parser.romNumSteps(),
				      obs, solver);
  }
  else if (parser.romNonlinearSolverType() == "GaussNewtonQR")
  {

    using mat_t = typename RomProblemType::traits::lspg_jacobian_type;
    using qr_solver_t = pressio::qr::QRSolver<mat_t, pressio::qr::Householder>;
    qr_solver_t qrSolver;
    auto solver = pnlins::create_gauss_newtonQR(problem, romState, qrSolver);
    solver.setStoppingCriterion(pnlins::Stop::WhenGradientAbsoluteNormBelowTolerance);

    const auto tol = parser.romNonlinearSolverTolerance();
    solver.setTolerance(tol);

    pode::advance_n_steps_and_observe(problem, romState, 0.,
				      parser.romTimeStepSize(),
				      parser.romNumSteps(),
				      obs, solver);
  }

  else if (parser.romNonlinearSolverType() == "LevenbergMarquardt")
  {
    using hessian_t = Eigen::Matrix<scalar_t, -1, -1>;
    using solver_tag = pressio::linearsolvers::direct::HouseholderQR;
    using linear_solver_t = pressio::linearsolvers::Solver<solver_tag, hessian_t>;
    linear_solver_t linSolverObj;

    auto solver = pnlins::create_levenberg_marquardt(problem, romState, linSolverObj);
    solver.setUpdatingCriterion(pnlins::Update::LMSchedule2);
    solver.setStoppingCriterion(pnlins::Stop::WhenGradientAbsoluteNormBelowTolerance);

    const auto tol = parser.romNonlinearSolverTolerance();
    solver.setTolerance(tol);

    pode::advance_n_steps_and_observe(problem, romState, 0.,
				      parser.romTimeStepSize(),
				      parser.romNumSteps(),
				      obs, solver);
  }
}

template<class AppObjType, class ParserType>
void run_lspg_hyperreduced(const AppObjType & appObj,
			   const ParserType & parser)
{
  // namespace pode = pressio::ode;
  // namespace prom = pressio::rom;
  // namespace pnlins = pressio::nonlinearsolvers;
  // namespace plspg  = pressio::rom::lspg;

  // using scalar_t = typename AppObjType::scalar_type;
  // using fom_state_t  = typename AppObjType::state_type;

  // // get initial condition
  // fom_state_t fomIC = appObj.initialCondition();

  // // podFile is file with the modes on the FULL mesh
  // const auto & podFile = parser.basisFile();
  // auto phiFull = create_colmajor_matrix_and_load_from_binary_with_extents<scalar_t>(podFile,
  // 										    parser.romSize());

  // // construct decoder on the stencil mesh
  // auto phiOnStencil = create_basis_on_stencil_mesh<scalar_t>(phiFull,
  // 							     parser.meshDir(),
  // 							     parser.numDofsPerCell());
  // auto decoder = prom::create_time_invariant_linear_decoder<fom_state_t>(phiOnStencil);

  // // create rom state
  // using rom_state_t = Eigen::Matrix<scalar_t, -1, 1>;
  // rom_state_t romState(parser.romSize());
  // romState = phiOnStencil.transpose() * fomIC;

  // // reference state
  // fom_state_t fomReferenceState(fomIC);
  // fomReferenceState.setZero();

  // const auto sampleMeshGidsFile  = parser.meshDir() + "/sample_mesh_gids.dat";
  // const auto stencilMeshGidsFile = parser.meshDir() + "/stencil_mesh_gids.dat";
  // auto sampleMeshGids  = read_gids_ascii_file(sampleMeshGidsFile);
  // auto stencilMeshGids = read_gids_ascii_file(stencilMeshGidsFile);

  // using combiner_t = HypRedUpdater<scalar_t>;
  // combiner_t combiner(sampleMeshGids, stencilMeshGids, parser.numDofsPerCell());

  // const auto odeScheme = parser.romOdeScheme();
  // auto problem = plspg::create_hyperreduced_unsteady_problem(odeScheme, appObj,
  // 							     decoder, romState,
  // 							     fomReferenceState,
  // 							     combiner);

  // pick_solver_and_run(parser, problem, romState);
}

template<class AppObjType, class ParserType>
void run_lspg_default(const AppObjType & appObj,
		      const ParserType & parser)
{
  // namespace prom = pressio::rom;
  // namespace plspg  = pressio::rom::lspg;

  // using scalar_t = typename AppObjType::scalar_type;
  // using fom_state_t  = typename AppObjType::state_type;

  // // get initial condition
  // fom_state_t fomIC = appObj.initialCondition();

  // // create decoder
  // const auto & podFile = parser.basisFile();
  // auto phi = create_colmajor_matrix_and_load_from_binary_with_extents<scalar_t>(podFile,
  // 										parser.romSize());
  // auto decoder = prom::create_time_invariant_linear_decoder<fom_state_t>(phi);

  // // reference state
  // fom_state_t fomReferenceState(fomIC);
  // //auto fomReferenceState = read_ref_state<scalar_t, fom_state_t>("ref_state.txt");
  // fomReferenceState.setZero();

  // // create rom state
  // using rom_state_t = Eigen::Matrix<scalar_t, -1, 1>;
  // rom_state_t romState(parser.romSize());
  // //romState = phi.transpose() *(fomIC - fomReferenceState);
  // romState = phi.transpose() * fomIC;

  // const auto odeScheme = parser.romOdeScheme();
  // auto problem = plspg::create_default_unsteady_problem(odeScheme, appObj,
  // 							decoder, romState,
  // 							fomReferenceState);

  // pick_solver_and_run(parser, problem, romState);
}

#endif
