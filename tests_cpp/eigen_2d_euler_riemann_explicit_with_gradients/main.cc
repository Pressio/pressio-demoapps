#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "../observer.hpp"

namespace pda = pressiodemoapps;

template<class MeshType, class T>
void write_gradients_fo_file(const MeshType & mesh,
			     const T & o,
			     const std::string & fileId)
{
  const std::string ssStr = std::to_string(mesh.stencilSize());
  const auto & G = mesh.graph();
  const auto & rowsCellsOnBD = mesh.graphRowsOfCellsStrictlyOnBd();
  std::ofstream file; file.open("grad_result_" + fileId + ".txt");

  auto toFile = [&](int cellGID, auto const & faceIn){
    file << faceIn.centerCoordinates[0] << " "
	 << faceIn.centerCoordinates[1] << " "
	 << std::setprecision(14)
	 << faceIn.normalGradient[0] << " "
	 << faceIn.normalGradient[1] << " "
	 << faceIn.normalGradient[2] << " "
	 << faceIn.normalGradient[3] << " "
	 << faceIn.normalDirection
	 << std::endl;
  };

  for (auto rowInd : rowsCellsOnBD){
    const bool bL = mesh.cellHasLeftFaceOnBoundary2d(rowInd);
    const bool bF = mesh.cellHasFrontFaceOnBoundary2d(rowInd);
    const bool bR = mesh.cellHasRightFaceOnBoundary2d(rowInd);
    const bool bB = mesh.cellHasBackFaceOnBoundary2d(rowInd);
    const int cellGID = G(rowInd, 0);

    if (bL){
      auto & face = o.queryFace(cellGID, pda::FacePosition::Left);
      toFile(cellGID, face);
    }
    if (bF){
      auto & face = o.queryFace(cellGID, pda::FacePosition::Front);
      toFile(cellGID, face);
    }
    if (bR){
      auto & face = o.queryFace(cellGID, pda::FacePosition::Right);
      toFile(cellGID, face);
    }
    if (bB){
      auto & face = o.queryFace(cellGID, pda::FacePosition::Back);
      toFile(cellGID, face);
    }
  }
  file.close();
}


template<class StateType>
void write_state_to_file(const StateType & state,
			 const std::string & filename)
{
  std::ofstream file1; file1.open(filename);
  for (int i=0;i <state.size(); ++i){ file1 << std::setprecision(14) << state(i) << '\n'; }
  file1.close();
}

int main()
{

  const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
#ifdef USE_WENO5
  const auto order   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
  const auto order   = pda::InviscidFluxReconstruction::Weno3;
#else
  const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

  const auto probId = pda::Euler2d::Riemann;
  auto appObj       = pda::create_problem_eigen(meshObj, probId, order, 2);

  constexpr int nDofPerCell = 4;
  pda::GradientEvaluator<decltype(meshObj), nDofPerCell> grads(meshObj);

  const auto dt = 0.01;
  const auto Nsteps = pressio::ode::StepCount(0.5/dt);
  Eigen::VectorXd state = appObj.initialCondition();

  write_state_to_file(state, "state_init.txt");
  grads(state, nDofPerCell);
  write_gradients_fo_file(meshObj, grads, "init");

  auto stepperObj = pressio::ode::create_ssprk3_stepper(appObj);
  pressio::ode::advance_n_steps(stepperObj, state, 0., dt, Nsteps);

  write_state_to_file(state, "state_final.txt");
  grads(state, nDofPerCell);
  write_gradients_fo_file(meshObj, grads, "final");

  return 0;
}
