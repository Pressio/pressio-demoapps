#ifndef PARSERS_HPP_
#define PARSERS_HPP_

#include "yaml-cpp/parser.h"
#include "yaml-cpp/yaml.h"

template<class T = void>
pressio::ode::StepScheme string_to_ode_scheme(const std::string & schemeString)
{

  if (schemeString == "ForwardEuler"){
    return pressio::ode::StepScheme::ForwardEuler;
  }
  else if (schemeString == "RungeKutta4"){
    return pressio::ode::StepScheme::RungeKutta4;
  }
  else if (schemeString == "SSPRungeKutta3"){
    return pressio::ode::StepScheme::SSPRungeKutta3;
  }
  else if (schemeString == "BDF1"){
    return pressio::ode::StepScheme::BDF1;
  }
  else if (schemeString == "CrankNicolson"){
    return pressio::ode::StepScheme::CrankNicolson;
  }
  else if (schemeString == "BDF2"){
    return pressio::ode::StepScheme::BDF2;
  }
  else{
    throw std::runtime_error("string_to_ode_scheme: Invalid odeScheme");
  }

}

template <typename ScalarType>
struct ParserCommon
{
private:
  std::string meshDirPath_ = "empty";
  int numDofsPerCell_ = {};

public:
  auto meshDir() const{ return meshDirPath_; }
  auto numDofsPerCell() const { return numDofsPerCell_; }

  ParserCommon() = delete;

  ParserCommon(YAML::Node & parentNode)
  {
    const auto node = parentNode["general"];
    if (node)
      {
	auto entry = "meshDir";
	if (node[entry]) meshDirPath_ = node[entry].as<std::string>();
	else throw std::runtime_error("Empty meshDir");

	entry = "numDofsPerCell";
	if (node[entry]) numDofsPerCell_ = node[entry].as<int>();
	else throw std::runtime_error("Empty numDofsPerCell");
      }
    else{
      throw std::runtime_error("General section in yaml input is mandatory!");
    }

    this->validate();
    this->describe();
  }

  void validate() const{
    if (meshDirPath_=="") {
      throw std::runtime_error("meshdir is empty string");
    }

    if (numDofsPerCell_ <= 0) {
      throw std::runtime_error("numDofsPerCell must be > 0");
    }
  }

  void describe() const{
    std::cout << "\nmeshDir = "	<< meshDirPath_	  << " \n"
	      << "numDofPerCell = " << numDofsPerCell_ << " \n";
  }
};

template <typename ScalarType>
struct ParserFomRomData
{
  bool isFom_                  = false;
  bool isRom_                  = false;
  ScalarType dt_	       = {};
  ScalarType finalTime_	       = {};
  int numSteps_                = {};
  int velocitySamplingFreq_    = 1;
  int stateSamplingFreq_       = {};
  std::string odeSchemeString_ = {};
  pressio::ode::StepScheme odeScheme_  = {};
  std::string nonlinSolverType_ = "void";
  ScalarType nonlinSolverTol_   = std::numeric_limits<ScalarType>::max();

  ParserFomRomData() = default;

  void parse(YAML::Node & node)
  {
    auto entry = "odeScheme";
    if (node[entry]) {
      odeSchemeString_ = node[entry].as<std::string>();
      odeScheme_ = string_to_ode_scheme(odeSchemeString_);
    }
    else{
      throw std::runtime_error("Input: missing odeScheme");
    }

    entry = "dt";
    if (node[entry]) dt_ = node[entry].as<ScalarType>();
    else throw std::runtime_error("Input: missing dt");

    entry = "finalTime";
    if (node[entry]) finalTime_ = node[entry].as<ScalarType>();
    else throw std::runtime_error("Input: missing final time");

    entry = "stateSamplingFreq";
    if (node[entry]) stateSamplingFreq_ = node[entry].as<int>();
    else throw std::runtime_error("Input: missing state sampling freq");

    entry = "velocitySamplingFreq";
    if (node[entry]) velocitySamplingFreq_ = node[entry].as<int>();

    entry = "nonlinearSolverType";
    if (node[entry]){
      nonlinSolverType_ = node[entry].as<std::string>();
    }
    else if (pressio::ode::is_implicit_scheme(odeScheme_)){
      // solvers are required only when doing implicit time integration
      throw std::runtime_error("Input: missing solver name");
    }

    entry = "nonlinearSolverTol";
    if (node[entry]){
      nonlinSolverTol_ = node[entry].as<ScalarType>();
    }
    else if (pressio::ode::is_implicit_scheme(odeScheme_)){
      // solvers are required only when doing implicit time integration
      throw std::runtime_error("Input: missing solver tolerance");
    }

    numSteps_ = static_cast<int>(finalTime_/dt_);
    this->validate();
  }

  void validate() const{
    if (dt_ <= static_cast<ScalarType>(0)) {
      throw std::runtime_error("Cannot have dt <= 0");
    }

    if (finalTime_ <= static_cast<ScalarType>(0)) {
      throw std::runtime_error("Cannot have finalT <= 0");
    }

    if (stateSamplingFreq_<=0) {
      throw std::runtime_error("Cannot have state sampling freq <= 0");
    }
  }

};


template <typename ScalarType>
struct ParserFom : private ParserFomRomData<ScalarType>
{
  using base_t = ParserFomRomData<ScalarType>;

public:
  auto doingFom()                    const { return base_t::isFom_; }
  auto fomTimeStepSize()	     const { return base_t::dt_; }
  auto fomNumSteps()		     const { return base_t::numSteps_; }
  auto fomVelocitySamplingFreq()     const { return base_t::velocitySamplingFreq_; }
  auto fomStateSamplingFreq()	     const { return base_t::stateSamplingFreq_; }
  auto fomOdeScheme()		     const { return base_t::odeScheme_; }
  auto fomNonlinearSolverType()	     const { return base_t::nonlinSolverType_; }
  auto fomNonlinearSolverTolerance() const { return base_t::nonlinSolverTol_; }

  ParserFom() = delete;
  ParserFom(YAML::Node & parentNode) : base_t()
  {
    auto fomNode = parentNode["fom"];
    if (fomNode){
      base_t::isFom_ = true;
      base_t::isRom_ = false;
      base_t::parse(fomNode);
      this->describe();
    }
  }

  void describe() const{
    std::cout << "\nfomOdeScheme = " << base_t::odeSchemeString_ << " \n"
	      << "fomTimeStep = "    << base_t::dt_ << " \n"
	      << "fomFinalTime = "   << base_t::finalTime_ << " \n"
	      << "fomNumSteps = "    << base_t::numSteps_ << " \n"
	      << "fomStateSampFreq = "    << base_t::stateSamplingFreq_ << " \n"
	      << "fomVelocitySampFreq = " << base_t::velocitySamplingFreq_ << " \n";

    if (pressio::ode::is_implicit_scheme(base_t::odeScheme_)){
      std::cout << "fomNonlinearSolverType = " << base_t::nonlinSolverType_ << " \n"
		<< "fomNonlinearSolverTol = "  << base_t::nonlinSolverTol_ << " \n";
    }
  }
};

template <typename ScalarType>
struct ParserRom : private ParserFomRomData<ScalarType>
{
  using base_t = ParserFomRomData<ScalarType>;
  std::string romKind_ = "empty";
  int romSize_= {};
  std::string fullMeshBasisFileName_ = "empty";
  std::string projectorFileName_ = "empty";
  std::string referenceStateFileName_ = "empty";
  std::string romInitialStateFileName_ = "empty";
  std::string sampleMeshGidsFileName_ = "empty";
  std::string stencilMeshGidsFileName_ = "empty";

public:
  auto doingRom()                    const { return base_t::isRom_; }
  auto romTimeStepSize()	     const { return base_t::dt_; }
  auto romNumSteps()		     const { return base_t::numSteps_; }
  auto romVelocitySamplingFreq()     const { return base_t::velocitySamplingFreq_; }
  auto romStateSamplingFreq()	     const { return base_t::stateSamplingFreq_; }
  auto romOdeScheme()		     const { return base_t::odeScheme_; }
  auto romNonlinearSolverType()	     const { return base_t::nonlinSolverType_; }
  auto romNonlinearSolverTolerance() const { return base_t::nonlinSolverTol_; }

  auto romKind()		     const{ return romKind_; }
  int romSize()			     const{ return romSize_; }
  auto fullMeshBasisFile()           const{ return fullMeshBasisFileName_; }
  auto projectorFile()               const{ return projectorFileName_; }
  auto referenceStateFile()	     const{ return referenceStateFileName_; }
  auto romInitialStateFile()	     const{ return romInitialStateFileName_; }
  auto sampleMeshGidsFileName()	     const{ return sampleMeshGidsFileName_; }
  auto stencilMeshGidsFileName()     const{ return stencilMeshGidsFileName_; }

  ParserRom() = delete;
  ParserRom(YAML::Node & parentNode) : base_t()
  {
    auto romNode = parentNode["rom"];
    if (romNode){
      base_t::isFom_ = false;
      base_t::isRom_ = true;
      base_t::parse(romNode);

      auto entry = "numModes";
      if (romNode[entry]) romSize_ = romNode[entry].as<int>();
      else throw std::runtime_error("Input: ROM: missing # of modes");

      entry = "podFile";
      if (romNode[entry]) fullMeshBasisFileName_ = romNode[entry].as<std::string>();
      else throw std::runtime_error("Input: ROM: missing podFile");

      entry = "algoName";
      if (romNode[entry]) romKind_ = romNode[entry].as<std::string>();
      else throw std::runtime_error("Input: ROM: missing algoName");

      if (!validAlgoNames(romKind_)){
	throw std::runtime_error("Input: ROM: invalid algoName");
      }

      entry = "referenceStateFile";
      if (romNode[entry]) referenceStateFileName_ = romNode[entry].as<std::string>();
      else throw std::runtime_error("Input: ROM: missing referenceStateFile");

      entry = "romInitialStateFile";
      if (romNode[entry]) romInitialStateFileName_ = romNode[entry].as<std::string>();
      else throw std::runtime_error("Input: ROM: missing romInitialStateFile");

      const auto projEntry = "projectorFile";
      const auto projFound = romNode[projEntry];
      if (romKind_ == "GalerkinCollocation" ||
	  romKind_ == "GalerkinGappy" ||
	  romKind_ == "GalerkinMaskedCollocation" ||
	  romKind_ == "GalerkinMaskedGappy")
	{
	  if (!projFound){
	    throw std::runtime_error("Input: ROM: missing projector");
	  }

	  projectorFileName_ = romNode[projEntry].as<std::string>();
	}

      const auto smGidsEntry  = "sampleMeshGidsFile";
      const auto stmGidsEntry = "stencilMeshGidsFile";
      const auto smGidsFound  = romNode[smGidsEntry];
      const auto stmGidsFound = romNode[stmGidsEntry];
      if (romKind_ == "GalerkinCollocation" ||
	  romKind_ == "GalerkinGappy" ||
	  romKind_ == "GalerkinMaskedCollocation" ||
	  romKind_ == "GalerkinMaskedGappy" ||
	  romKind_ == "LspgCollocation" ||
	  romKind_ == "LspgGappy")
	{
	  if (!smGidsFound || !stmGidsFound){
	    throw std::runtime_error("Input: ROM: missing sample and stencil gids files");
	  }

	  sampleMeshGidsFileName_  = romNode[smGidsEntry].as<std::string>();
	  stencilMeshGidsFileName_ = romNode[stmGidsEntry].as<std::string>();
	}

      this->describe();
    }
  }

  void describe() const{
    std::cout << "\nromKind_ = " << romKind_ << " \n"
	      << "romSize_ = "   << romSize_ << " \n"
	      << "fullMeshBasisFileName_ = " << fullMeshBasisFileName_ << " \n"
	      << "romOdeScheme = " << base_t::odeSchemeString_ << " \n"
	      << "romTimeStep = "  << base_t::dt_ << " \n"
	      << "romFinalTime = " << base_t::finalTime_ << " \n"
	      << "romNumSteps = "  << base_t::numSteps_ << " \n"
	      << "romStateSampFreq = " << base_t::stateSamplingFreq_ << " \n";

    if (romKind_ == "GalerkinCollocation" ||
	romKind_ == "GalerkinGappy" ||
	romKind_ == "GalerkinMaskedCollocation" ||
	romKind_ == "GalerkinMaskedGappy")
    {
      std::cout << "projectorFile = " << projectorFileName_ << " \n";
    }

    if (romKind_ == "GalerkinCollocation" ||
	romKind_ == "GalerkinGappy" ||
	romKind_ == "GalerkinMaskedCollocation" ||
	romKind_ == "GalerkinMaskedGappy" ||
	romKind_ == "LspgCollocation" ||
	romKind_ == "LspgGappy")
      {
	std::cout << "sampleMeshGidsFileName  = " << sampleMeshGidsFileName_ << " \n"
		  << "stencilMeshGidsFileName = " << stencilMeshGidsFileName_
		  << "\n";
      }

    if (pressio::ode::is_implicit_scheme(base_t::odeScheme_)){
      std::cout << "romNonlinearSolverType = " << base_t::nonlinSolverType_ << " \n"
		<< "romNonlinearSolverTol = "  << base_t::nonlinSolverTol_ << " \n";
    }
  }

private:
  bool validAlgoNames(const std::string & algoName) const{
    if (algoName != "GalerkinFull" &&
	algoName != "GalerkinCollocation" &&
	algoName != "GalerkinGappy" &&
	algoName != "GalerkinMaskedCollocation" &&
	algoName != "GalerkinMaskedGappy" &&
	algoName != "LspgFull" &&
	algoName != "LspgCollocation" &&
	algoName != "LspgGappy")
      {
	return false;
      }

    return true;
  }
};

template<class ScalarType>
class ParserTwoDimShallowWater
  : public ParserCommon<ScalarType>,
    public ParserFom<ScalarType>,
    public ParserRom<ScalarType>
{
  std::string inviscidFluxRec_ = "";
  ScalarType gravity_  = 9.8;
  ScalarType coriolis_ = -3.0;
  ScalarType pulseMag_ = 0.125;

public:
  ParserTwoDimShallowWater() = delete;

  ParserTwoDimShallowWater(YAML::Node & parentNode)
    : ParserCommon<ScalarType>(parentNode),
      ParserFom<ScalarType>(parentNode),
      ParserRom<ScalarType>(parentNode)
  {
    // check if we are doing FOM or ROM
    const auto isFom = this->doingFom();
    const auto isRom = this->doingRom();
    if (isFom && isRom){
      throw std::runtime_error("Seems like both FOM and ROM are on, this is not possible");
    }

    // depending on what is on, find reconstruction
    {
      const YAML::Node node = isFom ? parentNode["fom"] : parentNode["rom"] ;
      if (node["inviscidFluxReconstruction"]){
	inviscidFluxRec_ = node["inviscidFluxReconstruction"].as<std::string>();
      }
      else{
	throw std::runtime_error("Input: SWE: missing inviscid reconstruction");
      }
    }

    // find coefficients
    {
      const auto node = parentNode["physicalCoefficients"];
      if (node)
	{
	  auto entry = "gravity";
	  if (node[entry]) gravity_ = node[entry].as<ScalarType>();

	  entry = "coriolis";
	  if (node[entry]) coriolis_ = node[entry].as<ScalarType>();

	  entry = "pulsemag";
	  if (node[entry]) pulseMag_ = node[entry].as<ScalarType>();

	}
      else{
	throw std::runtime_error("SWE: missing coefficients section in yaml input");
      }
    }

    this->validate();
    this->describe();
  }

  auto inviscidFluxReconstruction() const{ return inviscidFluxRec_; }
  auto gravity() const{ return gravity_; }
  auto coriolis() const{ return coriolis_; }
  auto pulseMag() const{ return pulseMag_; }

private:
  void validate() const{
    if (inviscidFluxRec_ != "FirstOrder" &&
	inviscidFluxRec_ != "Weno3" &&
	inviscidFluxRec_ != "Weno5")
    {
      throw std::runtime_error
	("Input: SWE: inviscidFluxReconstruction must be: FirstOrder, Weno3, or Weno5");
    }

    if (gravity_<=0.) {
      throw std::runtime_error("Input: SWE: Cannot have negative gravity");
    }
    if (pulseMag_<=0.) {
      throw std::runtime_error("Input: SWE: Cannot have negative pulse magnitude");
    }
  }

  void describe() const{
    std::cout << "\ninviscidFluxReconstruction = " << inviscidFluxRec_ << " \n"
	      << "gravity = "   << gravity_ << " \n"
	      << "coriolis  = " << coriolis_ << " \n"
	      << "pulse_mag = " << pulseMag_ << " \n";
  }
};


template<class ScalarType>
class ParserTwoDimGrayScott
  : public ParserCommon<ScalarType>,
    public ParserFom<ScalarType>,
    public ParserRom<ScalarType>
{
  ScalarType diffusionA_ = {};
  ScalarType diffusionB_ = {};
  ScalarType feedRate_   = {};
  ScalarType killRate_   = {};

public:
  ParserTwoDimGrayScott() = delete;

  ParserTwoDimGrayScott(YAML::Node & parentNode)
    : ParserCommon<ScalarType>(parentNode),
      ParserFom<ScalarType>(parentNode),
      ParserRom<ScalarType>(parentNode)
  {
    // check if we are doing FOM or ROM
    const auto isFom = this->doingFom();
    const auto isRom = this->doingRom();
    if (isFom && isRom){
      throw std::runtime_error("Seems like both FOM and ROM are on, this is not possible");
    }

    // find coefficients
    {
      const auto node = parentNode["physicalCoefficients"];
      if (node)
	{
	  auto entry = "diffusionA";
	  if (node[entry]) diffusionA_ = node[entry].as<ScalarType>();
	  else throw std::runtime_error("missing diffusionA");

	  entry = "diffusionB";
	  if (node[entry]) diffusionB_ = node[entry].as<ScalarType>();
	  else throw std::runtime_error("missing diffusionB");

	  entry = "feedRate";
	  if (node[entry]) feedRate_ = node[entry].as<ScalarType>();
	  else throw std::runtime_error("missing feedRate");

	  entry = "killRate";
	  if (node[entry]) killRate_ = node[entry].as<ScalarType>();
	  else throw std::runtime_error("missing killRate");
	}
      else{
	throw std::runtime_error("2d GrayScott: missing coefficients section in yaml input");
      }
    }

    this->describe();
  }

  auto diffusionA() const{ return diffusionA_; }
  auto diffusionB() const{ return diffusionB_; }
  auto feedRate() const{ return feedRate_; }
  auto killRate() const{ return killRate_; }

private:
  void describe() const{
    std::cout << "\ndiffusionA = " << diffusionA_ << " \n"
	      << "diffusionB = " << diffusionB_ << " \n"
	      << "killRate = " << killRate_ << " \n"
	      << "feedRate = " << feedRate_ << " \n";
  }
};



// template<class ScalarType>
// class ParserTwoDimReacDiff
//   : public ParserCommon<ScalarType>,
//     public ParserRom<ScalarType>
// {
//   ScalarType diffusion_ = 0.01;
//   ScalarType reaction_  = 0.02;

// public:
//   ParserTwoDimReacDiff() = delete;

//   ParserTwoDimReacDiff(YAML::Node & node0)
//     : ParserCommon<ScalarType>(node0),
//       ParserRom<ScalarType>(node0)
//   {
//     const auto node = node0["parameters"];
//     if (node)
//     {
//       auto entry = "diffusion";
//       if (node[entry]) diffusion_ = node[entry].as<ScalarType>();
//       else throw std::runtime_error("missing diffusion");

//       entry = "reaction";
//       if (node[entry]) reaction_ = node[entry].as<ScalarType>();
//       else throw std::runtime_error("missing reaction");
//     }
//     else{
//       throw std::runtime_error("General section in yaml input is mandatory!");
//     }
//     this->validate();
//     this->describe();
//   }

//   auto diffusion() const{ return diffusion_; }
//   auto reaction() const{ return reaction_; }

// private:
//   void validate() const{
//     if (diffusion_<=0.) {
//       throw std::runtime_error("Cannot have negative diffusion");
//     }
//     if (reaction_<=0.){
//       throw std::runtime_error("Cannot have negative reaction");
//     }
//   }

//   void describe() const{
//     std::cout << "\ndiffusion = " << diffusion_ << " \n"
// 	      << "reaction = " << reaction_ << " \n";
//   }
// };


#endif
