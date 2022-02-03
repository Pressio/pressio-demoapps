
#ifndef RUN_GALERKIN_HPP_
#define RUN_GALERKIN_HPP_

#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressio/rom_galerkin.hpp"
#include "observer.hpp"
#include "rom_shared.hpp"

template<class ScalarType, class FomVeloType>
struct Masker
{
  std::vector<int> indices_;
  int numDofsPerCell_ = {};

  template<class ParserType>
  Masker(const ParserType & parser)
    : numDofsPerCell_(parser.numDofsPerCell())
  {
    fill_gids_stdvector_from_ascii(parser.sampleMeshGidsFileName(), indices_);
  }

  FomVeloType createApplyMaskResult(const FomVeloType & /*unmasked_object*/) const{
    return FomVeloType( numDofsPerCell_*indices_.size() );
  }

  template<class TimeType>
  void operator()(const FomVeloType & unmasked_object,
                  const TimeType time,
                  FomVeloType & result) const
  {
    for (int i=0; i<indices_.size(); ++i){
      for (int k=0; k<numDofsPerCell_; ++k){
	const int r = i*numDofsPerCell_ + k;
	const int g = indices_[i]*numDofsPerCell_ + k;
	result(r) = unmasked_object(g);
      }
    }
  }
};

template<class ScalarType>
struct ArbitraryProjector
{
  using operator_type = Eigen::Matrix<ScalarType, -1, -1, Eigen::ColMajor>;
  operator_type P_;

  ArbitraryProjector(const std::string & projFile,
		     const int numModes)
  {
    P_ = create_colmajor_matrix_and_load_from_binary_with_extents<
      ScalarType, operator_type>(projFile, numModes);
  }

  void operator()(const Eigen::Matrix<ScalarType, -1, 1> & operand,
		  const ScalarType time,
		  Eigen::Matrix<ScalarType, -1, 1> & result) const
  {
    result = P_.transpose() * operand;
  }
};

template<class AppObjType, class ParserType>
void run_galerkin_hyperreduced(const AppObjType & appObj,
			       const ParserType & parser)
{
  namespace pode = pressio::ode;
  namespace prom = pressio::rom;
  using scalar_t    = typename AppObjType::scalar_type;

  std::cout << "run_galerkin_hyperreduced \n";

  // podFile is file with the modes on the FULL mesh
  auto phiFull = create_colmajor_matrix_and_load_from_binary_with_extents<scalar_t>(parser.fullMeshBasisFile(),
										    parser.romSize());
  // decoder on the stencil mesh
  auto phiOnStencil = create_basis_on_stencil_mesh(phiFull, parser);
  using fom_state_t = typename AppObjType::state_type;
  auto decoder = prom::create_time_invariant_linear_decoder<fom_state_t>(phiOnStencil);

  // rom state
  auto romState = create_rom_state_and_load_from_ascii<scalar_t>(parser.romInitialStateFile());

  // observer
  StateObserver<decltype(romState)> observer("rom_snaps.bin", parser.romStateSamplingFreq());

  // projector
  ArbitraryProjector<scalar_t> P(parser.projectorFile(), romState.size());

  // reference state
  auto fomReferenceState = create_reference_state_on_stencil_mesh<scalar_t>(parser);

  const auto odeScheme = parser.romOdeScheme();
  if (pressio::ode::is_explicit_scheme(odeScheme))
  {
    auto problem = prom::galerkin::create_hyperreduced_explicit_problem(
      odeScheme, appObj, decoder, romState, fomReferenceState, P);

    pode::advance_n_steps_and_observe(problem, romState, 0.,
				      parser.romTimeStepSize(),
				      parser.romNumSteps(),
				      observer);
  }else{
    throw std::runtime_error("Implicit Galerkin not impl yet");
  }

}

template<class AppObjType, class ParserType>
void run_galerkin_masked(const AppObjType & appObj,
			 const ParserType & parser)
{

  namespace pode = pressio::ode;
  namespace prom = pressio::rom;
  using scalar_t    = typename AppObjType::scalar_type;

  std::cout << "run_galerkin_masked \n";

  // podFile is file with the modes on the FULL mesh
  auto phiFull = create_colmajor_matrix_and_load_from_binary_with_extents<scalar_t>(parser.fullMeshBasisFile(),
										    parser.romSize());
  using fom_state_t = typename AppObjType::state_type;
  auto decoder = prom::create_time_invariant_linear_decoder<fom_state_t>(phiFull);

  // rom state
  auto romState = create_rom_state_and_load_from_ascii<scalar_t>(parser.romInitialStateFile());

  // observer
  StateObserver<decltype(romState)> observer("rom_snaps.bin", parser.romStateSamplingFreq());

  // projector
  ArbitraryProjector<scalar_t> P(parser.projectorFile(), romState.size());

  // masker
  using fom_velo_t  = typename AppObjType::velocity_type;
  Masker<scalar_t, fom_velo_t> masker(parser);

  // reference state
  auto fomReferenceState = create_full_reference_state_and_load_from_ascii<scalar_t>(parser.referenceStateFile());

  const auto odeScheme = parser.romOdeScheme();
  if (pressio::ode::is_explicit_scheme(odeScheme))
  {
    auto problem = prom::galerkin::create_masked_explicit_problem(
     odeScheme, appObj, decoder, romState, fomReferenceState, P, masker);

    pode::advance_n_steps_and_observe(problem, romState, 0.,
				      parser.romTimeStepSize(),
				      parser.romNumSteps(),
				      observer);
  }else{
    throw std::runtime_error("Implicit Galerkin not impl yet");
  }

}

template<class AppObjType, class ParserType>
void run_galerkin_default(const AppObjType & appObj,
			  const ParserType & parser)
{
  namespace pode = pressio::ode;
  namespace prom = pressio::rom;
  using scalar_t = typename AppObjType::scalar_type;

  std::cout << "run_galerkin_default \n";

  // decoder
  auto phi = create_colmajor_matrix_and_load_from_binary_with_extents<scalar_t>(parser.fullMeshBasisFile(),
										parser.romSize());
  using fom_state_t  = typename AppObjType::state_type;
  auto decoder = prom::create_time_invariant_linear_decoder<fom_state_t>(phi);

  // rom state
  auto romState = create_rom_state_and_load_from_ascii<scalar_t>(parser.romInitialStateFile());

  // observer
  StateObserver<decltype(romState)> observer("rom_snaps.bin", parser.romStateSamplingFreq());

  // reference state
  auto fomReferenceState = create_full_reference_state_and_load_from_ascii<scalar_t>(parser.referenceStateFile());

  const auto odeScheme = parser.romOdeScheme();
  if (pressio::ode::is_explicit_scheme(odeScheme))
  {
    auto problem = prom::galerkin::create_default_explicit_problem(
     odeScheme, appObj, decoder, romState, fomReferenceState);

    pode::advance_n_steps_and_observe(problem, romState, 0.,
				      parser.romTimeStepSize(),
				      parser.romNumSteps(),
				      observer);
  }else{
    throw std::runtime_error("Implicit Galerkin not impl yet");
  }

}

#endif
