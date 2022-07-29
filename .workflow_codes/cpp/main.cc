
#include "pressiodemoapps/euler1d.hpp"
#include "pressiodemoapps/diffusion_reaction.hpp"
#include "pressiodemoapps/swe2d.hpp"
#include "run_fom.hpp"
#include "run_galerkin.hpp"
#include "run_lspg.hpp"
#include "parsers.hpp"
#include "source.hpp"
#include <chrono>

template<class AppObjType, class ParserType>
void dispatch(const AppObjType & appObj,
	      const ParserType & parser)
{

  if (parser.doingFom()){
    run_fom(appObj, parser);
  }

  else if (parser.doingRom() &&
	   parser.romKind()=="GalerkinFull")
  {
    run_galerkin_default(appObj, parser);
  }

  else if (parser.doingRom() &&
	   (parser.romKind()=="GalerkinCollocation" ||
	    parser.romKind()=="GalerkinGappy"))
  {
    run_galerkin_hyperreduced(appObj, parser);
  }

  else if (parser.doingRom() &&
	   (parser.romKind()=="GalerkinMaskedCollocation"))
  {
    run_galerkin_masked(appObj, parser);
  }

  else if (parser.doingRom() &&
	   (parser.romKind()=="GalerkinMaskedGappy"))
  {
    run_galerkin_masked(appObj, parser);
  }

  else if (parser.doingRom() &&
	   parser.romKind()=="LspgFull")
  {
    run_lspg_default(appObj, parser);
  }

  else if (parser.doingRom() &&
	   (parser.romKind()=="LspgCollocation" ||
	    parser.romKind()=="LspgGappy"))
  {
    run_lspg_hyperreduced(appObj, parser);
  }

  else{
    throw std::runtime_error("Invalid branch");
  }
}

std::string check_and_get_inputfile(int argc, char *argv[])
{
  if (argc != 2){
    throw std::runtime_error("Call as: ./exe <path-to-inputfile>");
  }
  const std::string inputFile = argv[1];
  std::cout << "Input file: " << inputFile << "\n";
  assert( file_exists(inputFile) );
  return inputFile;
}

int main(int argc, char *argv[])
{
  using scalar_t = double;

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::info});

  const auto inputFile = check_and_get_inputfile(argc, argv);
  auto node            = YAML::LoadFile(inputFile);
  const auto genNode   = node["general"];
  if (!genNode) throw std::runtime_error("Missing general section in yaml input!");

  std::string problem = {};
  if (genNode["problem"]){
    problem = genNode["problem"].as<std::string>();
  }else{
    throw std::runtime_error("Missing key=problem in yaml input");
  }

  // start timer
  auto t1 = std::chrono::high_resolution_clock::now();

  if (problem == "2d_swe")
  {
    namespace pda  = pressiodemoapps;
    ParserTwoDimShallowWater<scalar_t> parser(node);
    const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen<scalar_t>(parser.meshDir());
    const auto probId  = pda::Swe2d::SlipWall;
    const auto scheme  = (parser.inviscidFluxReconstruction()=="FirstOrder")
      ? pda::InviscidFluxReconstruction::FirstOrder :
         (parser.inviscidFluxReconstruction()=="Weno3")
            ? pda::InviscidFluxReconstruction::Weno3
               : pda::InviscidFluxReconstruction::Weno5;

    const auto gravity  = parser.gravity();
    const auto coriolis = parser.coriolis();
    const auto pulseMag = parser.pulseMag();
    auto appObj = pda::create_problem_eigen(meshObj, probId, scheme,
					    gravity, coriolis, pulseMag);
    dispatch(appObj, parser);
  }

  if (problem == "2d_gs")
  {
    namespace pda  = pressiodemoapps;
    ParserTwoDimGrayScott<scalar_t> parser(node);
    const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen<scalar_t>(parser.meshDir());
    const auto probId  = pda::DiffusionReaction2d::GrayScott;
    const auto scheme  = pda::ViscousFluxReconstruction::FirstOrder;

    const auto diffusionA = parser.diffusionA();
    const auto diffusionB = parser.diffusionB();
    const auto feedRate   = parser.feedRate();
    const auto killRate   = parser.killRate();
    auto appObj = pda::create_problem_eigen(meshObj, probId, scheme,
					    diffusionA, diffusionB,
					    feedRate, killRate);
    dispatch(appObj, parser);
  }

  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration< double > fs = t2 - t1;
  std::cout << "elapsed " << fs.count() << std::endl;

  pressio::log::finalize();
  return 0;
}




// if (problem == "2d_rd")
// {
//   namespace pda  = pressiodemoapps;
//   ParserTwoDimReacDiff<scalar_t> parser(node);
//   const auto meshObj   = pda::load_cellcentered_uniform_mesh_eigen<scalar_t>(parser.meshDir());
//   const auto diffusion = parser.diffusion();
//   const auto reaction  = parser.reaction();
//   const auto scheme    = pda::ViscousFluxReconstruction::FirstOrder;
//   const auto probId    = pda::DiffusionReaction2d::ProblemA;
//   auto appObj = pda::create_problem_eigen(meshObj, probId, scheme,
// 					    mySource, diffusion, reaction);
//   dispatch(appObj, parser);
// }
