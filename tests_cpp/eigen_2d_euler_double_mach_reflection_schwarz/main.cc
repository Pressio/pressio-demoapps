
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "../observer.hpp"

using namespace std;
namespace pda  = pressiodemoapps;
namespace plog = pressio::log;
namespace pode = pressio::ode;
namespace pls  = pressio::linearsolvers;
namespace pnls = pressio::nonlinearsolvers;

// TODO: Schwarz classes don't need access to appVec, should be fine with meshVec

// forward declarations
void calc_neighbor_dims(int, int, int, int, vector<vector<int>> &);
template<class app_t, class graph_t> vector<graph_t> calc_exch_graph(
  const int, const int, const int, const vector<app_t> &, const vector<vector<int>> &);
template<class app_t, class state_t, class graph_t> void broadcast_bcState(
  const int, const int, const int, vector<state_t> &, vector<state_t> &,
  const vector<app_t> &, const vector<vector<int>> &, const vector<graph_t> &);
template<class mesh_t, class graph_t> vector<graph_t> calc_ghost_graph(
  const int, const int, const vector<mesh_t> &, const vector<vector<int>> &);
template<class state_t> array<double, 2> calcConvergence(
  const state_t, const state_t, const int);

int main()
{

  plog::initialize(pressio::logto::terminal);
  plog::setVerbosity({pressio::log::level::debug});

  // ---- user inputs

  string meshRoot = "/home/crwentl/research/runs/pressio/double_mach/meshes/mesh_2x2";
  const auto order  = pda::InviscidFluxReconstruction::FirstOrder;
  const auto probId = pda::Euler2d::DoubleMachReflection;
  const auto scheme = pode::StepScheme::CrankNicolson;
  const int numSteps = 100;
  const int convergeStepMax = 10;
  vector<double> dt(4, 0.002);
  // vector<double> dt(4, 0.0025); // TODO: solver just gives up on domIdx > 0 for this time step
  double abs_err_tol = 1e-13;
  double rel_err_tol = 1e-13;
  // ---- end user inputs

  // load mesh info
  // TODO: move this into a function that just initializes meshVec
  int dim;
  int ndomX, ndomY, ndomZ;
  int overlap;
  pda::impl::read_domain_info(meshRoot, dim, ndomX, ndomY, ndomZ, overlap);
  const int ndomains = ndomX * ndomY;

  // dt setup
  if (dt.size() == 1) {
    dt.resize(ndomains, dt[0]);
  } else {
    if (dt.size() != ndomains) {
      cerr << "dt.size() must be 1 or ndomains, exiting" << endl;
      exit(-1);
    }
  }

  // controller step setup
  double dtMax = *max_element(dt.begin(), dt.end());
  vector<int> controlIters(ndomains);
  double integral;
  for (int domIdx = 0; domIdx < dt.size(); ++domIdx) {
    double niters = dtMax / dt[domIdx];
    if (round(niters) == niters) {
      controlIters[domIdx] = int(round(niters));
    } else {
      cerr << "dt of domain " << domIdx << " (" << dt[domIdx] << ") is not an integer divisor of maximum dt (" << dtMax << ")" << endl;
      exit(-1);
    }
  }

  // aliases
  using mesh_t  = decltype( pda::load_cellcentered_uniform_mesh_eigen(meshRoot + "/domain_0") );
  using app_t   = decltype( pda::create_problem_eigen(std::declval<mesh_t>(), probId, order) );
  using state_t = typename app_t::state_type;
  using jacob_t = typename app_t::jacobian_type;
  using lin_solver_t = pls::Solver<pls::iterative::Bicgstab, jacob_t>;
  using stepper_t = decltype( pode::create_implicit_stepper(scheme, std::declval<app_t &>()) );
  using nonlin_solver_t = decltype( pnls::create_newton_raphson( std::declval<stepper_t &>(), std::declval<lin_solver_t&>()) );
  using obs_t = FomObserver<state_t>;
  using graph_t = typename mesh_t::graph_t;

  // problem vectors initialization
  lin_solver_t            linSolverObj;
  vector<mesh_t>          meshVec(ndomains);
  vector<int>             ncellsVec(ndomains);
  vector<app_t>           appVec(ndomains);
  vector<state_t>         stateVec(ndomains);
  vector<vector<state_t>> stateHistVec(ndomains); // intra-controller step history, for interpolation
  vector<stepper_t>       stepperVec;
  vector<nonlin_solver_t> nonlinSolverVec;
  vector<obs_t>           obsVec(ndomains);
  for (int domIdx = 0; domIdx < ndomains; ++domIdx)
  {
    meshVec[domIdx]  = pda::load_cellcentered_uniform_mesh_eigen(meshRoot + "/domain_" + to_string(domIdx));
    ncellsVec[domIdx] = meshVec[domIdx].nx() * meshVec[domIdx].ny(); // TODO: generalize to 3D

    appVec[domIdx]   = pda::create_problem_eigen(meshVec[domIdx], probId, order);
    stateVec[domIdx] = appVec[domIdx].initialCondition();
    for (int histIdx = 0; histIdx < controlIters[domIdx] + 1; ++histIdx) {  // includes initial controller step solution
      stateHistVec[domIdx].emplace_back(appVec[domIdx].initialCondition());
    }

    stepperVec.emplace_back( pode::create_implicit_stepper(scheme, appVec[domIdx]) );
    nonlinSolverVec.emplace_back( pnls::create_newton_raphson(stepperVec.back(), linSolverObj) );
    nonlinSolverVec[domIdx].setTolerance(1e-5);

    obsVec[domIdx] = obs_t("doubleMach2d_solution_" + to_string(domIdx) + ".bin", 1);
    obsVec[domIdx](::pressio::ode::StepCount(0), 0.0, stateVec[domIdx]);
  }

  // ++++++++ BOUNDARY SETUP +++++++++++

  // set up Schwarz boundary graphs
  const auto stencilSize = pda::reconstructionTypeToStencilSize(order);
  const int numDofPerCell = 4;
  const int bcStencilSize = (stencilSize - 1) / 2;
  const int bcStencilDof = bcStencilSize * numDofPerCell;
  const int maxDomNeighbors = 4; // generalize to 2*ndim

  // determine neighboring domain IDs,
  vector<vector<int>> exchDomIdVec(ndomains, vector<int>(maxDomNeighbors, -1));
  calc_neighbor_dims(ndomX, ndomY, ndomZ, bcStencilDof, exchDomIdVec);

  // set up boundary broadcast patterns
  const auto exchGraphVec = calc_exch_graph<app_t, graph_t>(ndomains, overlap, bcStencilSize, appVec, exchDomIdVec);

  // create stateBcVec, sized accordingly
  vector<state_t> stateBcVec(ndomains);
  for (int domIdx = 0; domIdx < ndomains; ++domIdx) {
    // TODO: count numDofStencilBc AS NEEDED; lowers memory usage, but makes communication more difficult
    const auto meshObj = meshVec[domIdx];
    int numDofStencilBc = 2 * bcStencilDof * (meshObj.nx() + meshObj.ny() + meshObj.nz());
    pda::resize(stateBcVec[domIdx], numDofStencilBc);
    for (int dof = 0; dof < stateBcVec[domIdx].size(); ++dof) {
      stateBcVec[domIdx](dof) = 0.0;
    }
  }

  // first communication, set internal pointer
  for (int domIdx = 0; domIdx < ndomains; ++domIdx) {
    broadcast_bcState<app_t, state_t, graph_t>(domIdx, numDofPerCell, bcStencilSize, stateVec, stateBcVec, appVec, exchDomIdVec, exchGraphVec);
    appVec[domIdx].setStateBc(&stateBcVec[domIdx]);
  }

  // graph linking cell index to bcState start index
  auto ghostGraphVec = calc_ghost_graph<mesh_t, graph_t>(ndomains, numDofPerCell, meshVec, exchDomIdVec);
  for (int domIdx = 0; domIdx < ndomains; ++domIdx) {
    appVec[domIdx].setGraphBc(&ghostGraphVec[domIdx]);
  }

  // +++++++++ SOLVE ++++++++++++

  // controller outer loop
  double time = 0.0;
  double abs_err, rel_err;
  vector<array<double, 2>> convergeVals(ndomains);
  for (int outerStep = 1; outerStep <= numSteps; ++outerStep)
  {
    cout << "Step " << outerStep << endl;

    // store initial step for resetting
    for (int domIdx = 0; domIdx < ndomains; ++domIdx) {
      stateHistVec[domIdx][0] = stateVec[domIdx];
    }

    // convergence
    int convergeStep = 0;
    while (convergeStep < convergeStepMax) {

      cerr << "Schwarz iteration " << convergeStep + 1 << endl;

      // domain loop
      for (int domIdx = 0; domIdx < ndomains; ++domIdx) {

        // reset to beginning of controller time
        auto timeDom = time;
        auto stepDom = outerStep * controlIters[domIdx];

        const auto dtDom = dt[domIdx];
        const auto dtWrap = pode::StepSize<double>(dtDom);

        // controller inner loop
        for (int innerStep = 0; innerStep < controlIters[domIdx]; ++innerStep) {

          const auto startTimeWrap = pode::StepStartAt<double>(timeDom);
          const auto stepWrap = pode::StepCount(stepDom);

          stepperVec[domIdx](stateVec[domIdx], startTimeWrap, stepWrap, dtWrap, nonlinSolverVec[domIdx]);

          // for last iteration, compute convergence criteria
          // important to do this before saving history, as stateHistVec still has last convergence loop's state
          if (innerStep == (controlIters[domIdx] - 1)) {
            convergeVals[domIdx] = calcConvergence<state_t>(stateVec[domIdx], stateHistVec[domIdx].back(), ncellsVec[domIdx]);
          }

          // store intra-step history
          stateHistVec[domIdx][innerStep + 1] = stateVec[domIdx];

          // set (interpolated) boundary conditions

          // update local step and time
          stepDom++;
          timeDom += dtDom;

        }

        // broadcast boundary conditions
        broadcast_bcState<app_t, state_t, graph_t>(domIdx, numDofPerCell, bcStencilSize, stateVec, stateBcVec, appVec, exchDomIdVec, exchGraphVec);

      }

      // check convergence for all domains, break if conditions met
      abs_err = 0.0;
      rel_err = 0.0;
      for (int domIdx = 0; domIdx < ndomains; ++domIdx) {
        abs_err += convergeVals[domIdx][0];
        rel_err += convergeVals[domIdx][1];
      }
      abs_err /= ndomains;
      rel_err /= ndomains;
      cerr << "Average abs err: " << abs_err << endl;
      cerr << "Average rel err: " << rel_err << endl;
      if ((rel_err < rel_err_tol) || (abs_err < abs_err_tol)) {
        break;
      }

      convergeStep++;

      // reset interior state if not converged
      for (int domIdx = 0; domIdx < ndomains; ++domIdx) {
        stateVec[domIdx] = stateHistVec[domIdx][0];
      }

    }

    // output and updates
    const auto stepWrap = pode::StepCount(outerStep);
    for (int domIdx = 0; domIdx < ndomains; ++domIdx) {
      obsVec[domIdx](stepWrap, time + dtMax, stateVec[domIdx]);
    }

    time += dtMax;

  }
  return 0;
}

void calc_neighbor_dims(
  int ndomX,
  int ndomY,
  int ndomZ,
  int bcStencilDof,
  vector<vector<int>> & exchDomIdVec)
{
  // TODO: extend to 3D

  const int ndomains = ndomX * ndomY * ndomZ;
  vector<int> numDofStencilBcVec(ndomains);

  for (int domIdx = 0; domIdx < ndomains; ++domIdx) {

    int i = domIdx % ndomX;
    int j = domIdx / ndomX;
    int neighborId;

    // left boundary
    if (i != 0) {
      neighborId = domIdx - 1;
      exchDomIdVec[domIdx][0] = neighborId;
    }

    // "front" boundary
    if (j != (ndomY - 1)) {
      neighborId = domIdx + ndomX;
      exchDomIdVec[domIdx][1] = neighborId;
    }

    // right boundary
    if (i != (ndomX - 1)) {
      neighborId = domIdx + 1;
      exchDomIdVec[domIdx][2] = neighborId;
    }

    // "back" boundary
    if (j != 0) {
      neighborId = domIdx - ndomX;
      exchDomIdVec[domIdx][3] = neighborId;
    }

    // TODO: for 3D, need "bottom" and "top"

  }
}

template<class app_t, class graph_t>
vector<graph_t> calc_exch_graph(
  const int ndomains,
  const int overlap,
  const int bcStencilSize,
  const vector<app_t> & appVec,
  const vector<vector<int>> & exchDomIdVec)
{
  // TODO: extend to 3D

  // BC cell indexing example
  // L/R is from bottom to top, F/B is from left to right
  // Trying to mix cell ordering and face ordering
  //                  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
  //                 ¦_(4, 1)_¦_(5, 1)_¦_(6, 1)_¦_(7, 1)_¦_(8, 1)_¦
  //  _ _ _ _ _ _ _ _¦_(4, 0)_¦_(5, 0)_¦_(6, 0)_¦_(7, 0)_¦_(8, 0)_¦_ _ _ _ _ _ _ _ _
  // ¦_(3,1)_¦_(3,0)_|________|________|________|________|________|_(12,1)_¦_(12,0)_|
  // ¦_(2,1)_¦_(2,0)_|________|________|________|________|________|_(11,1)_¦_(11,0)_|
  // ¦_(1,1)_¦_(1,0)_|________|________|________|________|________|_(10,1)_¦_(10,0)_|
  // ¦_(0,1)_¦_(0,0)_|________|________|________|________|________|_(9, 1)_¦_(9, 0)_|
  //                 ¦_(13,0)_¦_(14,0)_¦_(15,0)_¦_(16,0)_¦_(17,0)_¦
  //                 ¦_(13,1)_¦_(14,1)_¦_(15,1)_¦_(16,1)_¦_(17,1)_¦

  vector<graph_t> exchGraphVec(ndomains);

  for (int domIdx = 0; domIdx < ndomains; ++domIdx) {

    // this domain's mesh and dimensions
    const auto domMesh = appVec[domIdx].getMesh();
    const auto domGraph = domMesh.graph();
    int nx = domMesh.nx();
    int ny = domMesh.ny();

    const auto & x = domMesh.viewX();
    const auto & y = domMesh.viewY();

    // TODO: generalize to 3D
    pda::resize(exchGraphVec[domIdx], 2*nx + 2*ny, bcStencilSize);
    exchGraphVec[domIdx].setOnes();
    exchGraphVec[domIdx] *= -1;

    // loop through neighboring domains
    for (int neighIdx = 0; neighIdx < exchDomIdVec[domIdx].size(); ++neighIdx) {

      int neighDomIdx = exchDomIdVec[domIdx][neighIdx];
      if (neighDomIdx == -1) {
        continue;  // not a Schwarz BC
      }

      // neighboring domain mesh and dimensions
      const auto neighMesh = appVec[neighDomIdx].getMesh();
      int nxNeigh = neighMesh.nx();
      int nyNeigh = neighMesh.ny();

      int exchCellIdx;

      // east-west neighbors will share a row index
      // left
      if (neighIdx == 0) {
        if (ny != nyNeigh) {
          cerr << "Mesh y-dimension mismatch for domains " << domIdx << " v " << neighDomIdx << ": " << ny << " != " << nyNeigh << endl;
          exit(-1);
        }

        int bcCellIdx = 0; // left boundary is the start
        for (int yIdx = 0; yIdx < ny; ++yIdx) {
          for (int stencilIdx = 0; stencilIdx < bcStencilSize; ++stencilIdx) {
            // exchCellIdx = (nxNeigh * (yIdx + 1) - 1) - overlap - stencilIdx;
            exchCellIdx = (nxNeigh * (yIdx + 1) - 1) - overlap + 1 - stencilIdx; // toward boundary
            exchGraphVec[domIdx](bcCellIdx, stencilIdx) = exchCellIdx;
          }
          bcCellIdx++;
        }
      }

      // north-south neighbors will share a column index
      // "front"
      if (neighIdx == 1) {
        if (nx != nxNeigh) {
          cerr << "Mesh x-dimension mismatch for domains " << domIdx << " v " << neighDomIdx << ": " << nx << " != " << nxNeigh << endl;
          exit(-1);
        }

        int bcCellIdx = ny * bcStencilSize;  // skip left boundary indices
        for (int xIdx = 0; xIdx < nx; ++xIdx) {
          for (int stencilIdx = 0; stencilIdx < bcStencilSize; ++stencilIdx) {
            exchCellIdx = (overlap + stencilIdx) * nxNeigh + xIdx;
            exchGraphVec[domIdx](bcCellIdx, stencilIdx) = exchCellIdx;
          }
          bcCellIdx++;
        }
      }

      // right
      if (neighIdx == 2) {
        if (ny != nyNeigh) {
          cerr << "Mesh y-dimension mismatch for domains " << domIdx << " v " << neighDomIdx << ": " << ny << " != " << nyNeigh << endl;
          exit(-1);
        }

        int bcCellIdx = (ny + nx) * bcStencilSize; // skip left and "front" boundary indices
        for (int yIdx = 0; yIdx < ny; ++yIdx) {
          for (int stencilIdx = 0; stencilIdx < bcStencilSize; ++stencilIdx) {
            exchCellIdx = nxNeigh * yIdx + overlap + stencilIdx;
            // TODO: check for a different problem
            // exchCellIdx = nxNeigh * yIdx + overlap - 1 + stencilIdx; // toward boundary
            // exchCellIdx = nxNeigh * yIdx + overlap + 1 + stencilIdx; // away from boundary
            exchGraphVec[domIdx](bcCellIdx, stencilIdx) = exchCellIdx;
          }
          bcCellIdx++;
        }
      }

      // "back"
      if (neighIdx == 3) {
        if (nx != nxNeigh) {
          cerr << "Mesh x-dimension mismatch for domains " << domIdx << " v " << neighDomIdx << ": " << nx << " != " << nxNeigh << endl;
          exit(-1);
        }

        int bcCellIdx = (2*ny + nx) * bcStencilSize;  // skip left, "front", and right boundary indices
        for (int xIdx = 0; xIdx < nx; ++xIdx) {
          for (int stencilIdx = 0; stencilIdx < bcStencilSize; ++stencilIdx) {
            exchCellIdx = (nyNeigh - 1 - overlap - stencilIdx) * nxNeigh + xIdx;
            exchGraphVec[domIdx](bcCellIdx, stencilIdx) = exchCellIdx;
          }
          bcCellIdx++;
        }
      }

      // TODO: generalize to 3D

    } // neighbor loop
  } // domain loop

  return exchGraphVec;

}

template<class graph_t, class state_t>
void comm_stateBc(
  const int startIdx,
  const int endIdx,
  const int bcStencilSize,
  const int dofPerCell,
  const graph_t & exchGraph,
  state_t & bcState,
  const state_t & intState)
{

  int exchCellIdx;
  for (int bcCellIdx = startIdx; bcCellIdx < endIdx; ++bcCellIdx) {
    for (int stencilIdx = 0; stencilIdx < bcStencilSize; ++stencilIdx) {
      exchCellIdx = exchGraph(bcCellIdx, stencilIdx);
      for (int dof = 0; dof < dofPerCell; ++dof) {
        bcState((bcCellIdx + stencilIdx) * dofPerCell + dof) = intState(exchCellIdx * dofPerCell + dof);
      }
    }
  }

}

template<class app_t, class state_t, class graph_t>
void broadcast_bcState(
  const int domIdx,
  const int dofPerCell,
  const int bcStencilSize,
  vector<state_t> & stateVec,
  vector<state_t> & stateBcVec,
  const vector<app_t> & appVec,
  const vector<vector<int>> & exchDomIdVec,
  const vector<graph_t> & exchGraphVec)
{

  const auto domState = stateVec[domIdx];

  int exchCellIdx;
  int startIdx;
  int endIdx;
  for (int neighIdx = 0; neighIdx < exchDomIdVec[domIdx].size(); ++neighIdx) {

    int neighDomIdx = exchDomIdVec[domIdx][neighIdx];
    if (neighDomIdx == -1) {
      continue;  // not a Schwarz BC
    }

    const auto neighMesh = appVec[neighDomIdx].getMesh();
    int nxNeigh = neighMesh.nx();
    int nyNeigh = neighMesh.ny();
    auto * neighStateBc = &stateBcVec[neighDomIdx];
    const auto neighExchGraph = exchGraphVec[neighDomIdx];

    // this domain is the neighboring domain's left neighbor
    if (neighIdx == 2) {
      startIdx = 0;
      endIdx = nyNeigh;
      comm_stateBc(startIdx, endIdx, bcStencilSize, dofPerCell, neighExchGraph, *neighStateBc, domState);
    }

    // this domain is the neighboring domain's front neighbor
    if (neighIdx == 3) {
      startIdx = nyNeigh;
      endIdx = nyNeigh + nxNeigh;
      comm_stateBc(startIdx, endIdx, bcStencilSize, dofPerCell, neighExchGraph, *neighStateBc, domState);
    }

    // this domain is the neighboring domain's right neighbor
    if (neighIdx == 0) {
      startIdx = nyNeigh + nxNeigh;
      endIdx = 2 * nyNeigh + nxNeigh;
      comm_stateBc(startIdx, endIdx, bcStencilSize, dofPerCell, neighExchGraph, *neighStateBc, domState);
    }

    // this domain is the neighboring domain's back neighbor
    if (neighIdx == 1) {
      startIdx = 2 * nyNeigh + nxNeigh;
      endIdx = 2 * nyNeigh + 2 * nxNeigh;
      comm_stateBc(startIdx, endIdx, bcStencilSize, dofPerCell, neighExchGraph, *neighStateBc, domState);
    }

  }

}

template<class mesh_t, class graph_t>
vector<graph_t> calc_ghost_graph(
  const int ndomains,
  const int dofPerCell,
  const vector<mesh_t> & meshVec,
  const vector<vector<int>> & exchDomIdVec)
{

  vector<graph_t> ghostGraphVec(ndomains);

  for (int domIdx = 0; domIdx < ndomains; ++domIdx) {

    const auto meshObj = meshVec[domIdx];
    const auto intGraph = meshObj.graph();
    int nx = meshObj.nx();
    int ny = meshObj.ny();

    const auto & rowsBd = meshObj.graphRowsOfCellsNearBd();
    pda::resize(ghostGraphVec[domIdx], int(rowsBd.size()), 4); // TODO: generalize to 2 * ndim
    ghostGraphVec[domIdx].setOnes();
    ghostGraphVec[domIdx] *= -1;

    for (decltype(rowsBd.size()) it = 0; it < rowsBd.size(); ++it) {

      // ASK FRANCESCO: is there an instance when rowsBd[it] != intGraph(rowsBd[it], 0)?
      //  The indices appear to be identical
      // TODO: this is all totally wrong for higher order

      const auto smPt = rowsBd[it];
      const auto left0  = intGraph(smPt, 1);
      const auto front0 = intGraph(smPt, 2);
      const auto right0 = intGraph(smPt, 3);
      const auto back0  = intGraph(smPt, 4);

      int stencilIdx = 0; // first order
      int rowIdx = smPt / nx;
      int colIdx = smPt % nx;
      int bcCellIdx;

      if (left0 == -1) {
        if (exchDomIdVec[domIdx][0] != -1) {
          bcCellIdx = rowIdx;
          ghostGraphVec[domIdx](it, 0) = bcCellIdx * dofPerCell;
        }
      }

      if (front0 == -1) {
        if (exchDomIdVec[domIdx][1] != -1) {
          bcCellIdx = ny + colIdx;
          ghostGraphVec[domIdx](it, 1) = bcCellIdx * dofPerCell;
        }
      }

      if (right0 == -1) {
        if (exchDomIdVec[domIdx][2] != -1) {
          bcCellIdx = ny + nx + rowIdx;
          ghostGraphVec[domIdx](it, 2) = bcCellIdx * dofPerCell;
        }
      }

      if (back0 == -1) {
        if (exchDomIdVec[domIdx][3] != -1) {
          bcCellIdx = 2 * ny + nx + colIdx;
          ghostGraphVec[domIdx](it, 3) = bcCellIdx * dofPerCell;
        }
      }
      // TODO: extend to higher order, 3D

    } // boundary cell loop
  } // domain loop

  return ghostGraphVec;

}


template<class state_t>
array<double, 2> calcConvergence(
  const state_t state1,
  const state_t state2,
  const int numElems)
{
  // TODO: assumed to be an Eigen state, not sure how to generalize
  // TODO: compute convergence for each variable separately

  int numDOF = state1.size();
  if (state2.size() != numDOF) {
    cerr << "state1 size does not match state2 size, " << numDOF << " vs. " << state2.size() << endl;
    exit(-1);
  }
  double numDOFPerCellD = double(numDOF) / numElems;
  if (round(numDOFPerCellD) != numDOFPerCellD) {
    cerr << "Computed a non-integer number of variables: " << numDOFPerCellD << endl;
    exit(-1);
  }
  int numDOFPerCell = round(numDOFPerCellD);

  // absolute error
  double abs_err = (state1 - state2).squaredNorm();

  // handle edge cases for relative error
  double rel_err;
  double basenorm = state1.squaredNorm();
  if (basenorm > 0) {
    rel_err = abs_err / basenorm;
  }
  else {
    if (abs_err > 0) {
        rel_err = 1.0;
    }
    else {
        rel_err = 0.0;
    }
  }

  array<double, 2> errArr = {abs_err, rel_err};
  return errArr;

}