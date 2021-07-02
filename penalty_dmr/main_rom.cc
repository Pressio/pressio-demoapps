
#include "pressio_rom.hpp"
#include "euler2d.hpp"
#include "observer.hpp"
#include "CLI11.hpp"

using scalar_t = double;
using rom_state_t    = pressio::containers::Vector<Eigen::VectorXd>;
using rom_residual_t = pressio::containers::Vector<Eigen::VectorXd>;
using rom_jacobian_t = pressio::containers::DenseMatrix<Eigen::MatrixXd>;
using decoder_jac_t  = pressio::containers::MultiVector<Eigen::MatrixXd>;

template<class decoder_t, class fom_problem_t>
struct MyPenaltySystem
{
  using eigen_spmat_t = Eigen::SparseMatrix<scalar_t>;
  using fom_state_t   = typename decoder_t::fom_state_type;
  using penalty_t     = typename fom_problem_t::penalty_vector_t;
  using penalty_jac_t = typename fom_problem_t::penalty_jacobian_t;

  using scalar_type   = scalar_t;
  using state_type    = rom_state_t;
  using residual_type = state_type;
  using jacobian_type = rom_jacobian_t;

  MyPenaltySystem(const scalar_t alphaIn,
		  const fom_problem_t & fomObjectIn,
		  const int numModesIn,
		  const decoder_t & decoderIn,
		  const fom_state_t & fomReferenceStateIn,
		  const state_type & romStateAt_nIn)
    : m_alpha(alphaIn),
      m_numModes(numModesIn),
      m_numResidualCells(fomObjectIn.totalResidualCells()),
      m_numDofsTotal(fomObjectIn.totalDofStencilMesh()),
      m_fomObject(fomObjectIn),
      m_decoder(decoderIn),
      m_fomReferenceState(fomReferenceStateIn),
      m_romStateAt_n(romStateAt_nIn),
      m_fomStateAux(m_numDofsTotal),
      m_g(fomObjectIn.createPenalty()),
      m_gJac(fomObjectIn.createPenaltyJacobian())
  {}

  void setCurrentTime(const scalar_t timeIn){
    m_currentTime = timeIn;
  }

  void setTimeStepSize(const scalar_t dtIn){
    m_dt = dtIn;
  }

  residual_type createResidual() const {
    return residual_type(m_numModes);
  }

  jacobian_type createJacobian() const {
    return jacobian_type(m_numModes, m_numModes);
  }

  void residual(const state_type & romState,
                residual_type & residual) const
  {
    std::cout << "sys: residual: romStateExtent = " << romState.extent(0) << std::endl;
    std::cout << "sys: residual: residualExtent = " << residual.extent(0) << std::endl;

    using cnst = pressio::utils::constants<scalar_t>;
    constexpr auto one    = cnst::one();
    constexpr auto negOne = one*-1.;
    constexpr auto zero   = cnst::zero();

    const auto & phi = m_decoder.jacobianCRef();
    std::cout << "sys:res: phi size = " << phi.extent(0) << " " << phi.extent(1) << std::endl;

    // -----------------------
    // first part of residual
    // -----------------------
    // yFom_n = yRef + phi*x_n
    reconstructFomState(m_romStateAt_n);

    // compute f
    fom_state_t velo(m_fomObject.createVelocity());
    m_fomObject.velocity(*m_fomStateAux.data(), m_currentTime, *velo.data());

    // compute residual = phi^T f
    pressio::ops::product(pressio::transpose(), one,
			  phi, velo, zero, residual);

    // update residual to compute: (1/dt)*x_n+1 - (1/dt)*x_n - phi^T f(x_n)
    const auto dtInv = one/m_dt;
    pressio::ops::update(residual, negOne,
			 romState, dtInv,
			 m_romStateAt_n, -dtInv);

    // -1/dt * r
    pressio::ops::scale(residual, -dtInv);

    // -----------------------
    // second part of residual
    // -----------------------
    // yFom_n = yRef + phi*x_n+1
    reconstructFomState(romState);

    // compute g() and dg/dx
    m_fomObject.penalty(*m_fomStateAux.data(), m_currentTime, m_g);
    //std::cout << "sys:res: g size = " << m_g.size() << std::endl;

    m_fomObject.penaltyJacobian(*m_fomStateAux.data(), m_currentTime, m_gJac);
    // std::cout << "sys:res: m_gJac size = "
    // 	      << m_gJac.rows() << " " << m_gJac.cols() << std::endl;

    auto gJTg = m_gJac.transpose() * m_g;
    // std::cout << "sys:res: gJTg size = "
    // 	      << gJTg.rows() << " " << gJTg.cols() << std::endl;

    auto resPart2Native = -m_alpha * (phi.data()->transpose() * gJTg);
    residual_type resP2(resPart2Native);
    ::pressio::ops::update(residual, one, resP2, one);
  }

  void jacobian(const state_type & romState,
		jacobian_type& jacobian) const
  {
    std::cout << "sys: jacobian: romStateExtent = "
	      << romState.extent(0) << std::endl;
    std::cout << "sys: jacobian: jaobian = "
	      << jacobian.data()->rows() << " "
	      << jacobian.data()->cols() << std::endl;

    const auto & phi = m_decoder.jacobianCRef();
    const auto & nativePhi = *phi.data();
    auto & nativeJ = *jacobian.data();
    auto A = nativePhi.transpose() * nativePhi;
    nativeJ = A.transpose() * A;

    pressio::ops::scale(jacobian, m_alpha);
    auto diag = pressio::containers::diag(jacobian);
    const auto dtInvSq = 1./(m_dt*m_dt);
    for (int i=0; i<diag.extent(0); ++i){
      diag(i) += dtInvSq;
    }
    std::cout << *jacobian.data() << std::endl;
  }

private:
  void reconstructFomState(const state_type & romState) const
  {
    m_decoder.applyMapping(romState, m_fomStateAux);
    pressio::ops::update(m_fomStateAux, 1.0, m_fomReferenceState, 1.0);
  }

private:
  scalar_t m_currentTime = {};
  scalar_t m_dt = {};
  const scalar_t m_alpha = {};
  const int m_numModes = {};
  const int m_numResidualCells = {};
  const int m_numDofsTotal = {};

  const fom_problem_t & m_fomObject;
  const decoder_t & m_decoder;
  const fom_state_t & m_fomReferenceState;
  const state_type & m_romStateAt_n;
  mutable fom_state_t m_fomStateAux;
  mutable penalty_t     m_g;
  mutable penalty_jac_t m_gJac;

};

template <typename T>
auto readBasis(std::string filename, int romSize, T totDofStencilMesh)
{
  std::vector<std::vector<scalar_t>> A0;
  ::pressio::utils::readAsciiMatrixStdVecVec(filename, A0, romSize);

  decoder_jac_t phi(totDofStencilMesh, romSize);
  for (T i=0; i<totDofStencilMesh; i++){
    for (int j=0; j<romSize; j++){
      phi(i,j) = A0[i][j];
    }
  }
  return phi;
}

template<class fom_state_t>
auto createLinearDecoder(int numModes,
			 int totDofStencilMesh,
			 const std::string & podModesFilePath)
{
  using decoder_t = pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;
  auto phi = readBasis(podModesFilePath, numModes, totDofStencilMesh);
  return decoder_t(std::move(phi));
}

int main(int argc, char *argv[])
{
  namespace pda = pressiodemoapps;
  pressio::log::initialize(pressio::logto::fileAndTerminal, "log.txt");
  pressio::log::setVerbosity({pressio::log::level::trace, pressio::log::level::trace});

  {
    CLI::App app{"2D Double Mach reflection - ROM"};
    scalar_t dt	       = 0.0001;
    int sampleEvery    = 50;
    scalar_t finalTime = 0.25;
    int numModes       = 5;
    scalar_t alpha     = 0.5;
    std::string meshPath = ".";
    std::string podPath = ".";

    app.add_option("--podPath",      podPath,	  "Full path to pod modes");
    app.add_option("--meshPath",     meshPath,	  "Full path to mesh");
    app.add_option("--sampleFreq",   sampleEvery, "Frequency of sampling state");
    app.add_option("--dt",	     dt,	  "Default = 0.0001");
    app.add_option("-p,--numModes",  numModes,	  "Default = 5");
    app.add_option("-T,--finalTime", finalTime,   "Default = 0.25");
    app.add_option("-a,--alpha",     alpha,       "Default = 0.5");
    CLI11_PARSE(app, argc, argv);

    // create the FOM physics object
    const auto meshObj	       = pda::loadCellCenterUniformMeshEigen(meshPath);
    constexpr auto inviscidRec = pda::InviscidFluxReconstruction::Weno3;
    constexpr auto probId      = pda::Euler2d::DoubleMachReflection;
    auto fomObj		       = pda::createProblemEigen(meshObj, probId, inviscidRec);
    using fom_t		       = decltype(fomObj);

    // create decoder
    using fom_native_state_t = typename fom_t::state_type;
    using fom_state_t	     = pressio::containers::Vector<fom_native_state_t>;
    auto decoder = createLinearDecoder<fom_state_t>(numModes,
						    fomObj.totalDofStencilMesh(),
						    podPath+"/pod_modes.txt");

    // reference state (same as initial condition)
    fom_state_t fomReferenceState(fomObj.initialCondition());

    // rom states
    rom_state_t romStateAt_np1(numModes); // hat(x)_n+1
    rom_state_t romStateAt_n(numModes);   // hat(x)_n
    pressio::ops::set_zero(romStateAt_np1);
    pressio::ops::set_zero(romStateAt_n);

    // system object to provide to nonlinear solver
    using decoder_t = decltype(decoder);
    using system_type = MyPenaltySystem<decoder_t, fom_t>;
    system_type MySystemObj(alpha, fomObj, numModes, decoder,
			    fomReferenceState, romStateAt_n);

    // create solver
    using lin_solver_mat_t = typename system_type::jacobian_type;
    using lin_tag = pressio::solvers::linear::iterative::LSCG;
    using lin_solver_t   = pressio::solvers::linear::Solver<lin_tag, lin_solver_mat_t>;
    lin_solver_t linearSolverObj;
    auto NonLinSolver = pressio::solvers::nonlinear::createNewtonRaphson(MySystemObj,
									 romStateAt_np1,
									 linearSolverObj);

    constexpr auto stopwhen = pressio::solvers::nonlinear::stop::afterMaxIters;
    NonLinSolver.setStoppingCriterion(stopwhen);
    NonLinSolver.setMaxIterations(1);

    // loop over time
    RomObserver<rom_state_t> RomObs("rom_dmr_solution.bin", sampleEvery);
    const auto Nsteps = finalTime/dt;
    auto time = 0.0;
    for (int step=1; step<=Nsteps; ++step)
      {
	MySystemObj.setCurrentTime(time);
	MySystemObj.setTimeStepSize(dt);

	NonLinSolver.solve(MySystemObj, romStateAt_np1);
	time = static_cast<double>(step) * dt;
	RomObs(step, time, romStateAt_n);

	pressio::ops::deep_copy(romStateAt_n, romStateAt_np1);
      }

  }

  return 0;
}
