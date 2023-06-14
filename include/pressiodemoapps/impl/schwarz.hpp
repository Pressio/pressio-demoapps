

#ifndef PRESSIODEMOAPPS_SCHWARZ_HPP_
#define PRESSIODEMOAPPS_SCHWARZ_HPP_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "pressio/ode_steppers_implicit.hpp"

using namespace std;

namespace pressiodemoapps{ namespace impl {

  namespace pda = pressiodemoapps;
  namespace pode = pressio::ode;
  namespace pls  = pressio::linearsolvers;
  namespace pnls = pressio::nonlinearsolvers;

  template<
    class prob_t,
    class order_t,
    class scheme_t
  >
  class SchwarzDecomp
  {

    public:

      using mesh_t         = decltype( pda::load_cellcentered_uniform_mesh_eigen(declval<string&>()) );
      using graph_t        = typename mesh_t::graph_t;
      using app_t          = decltype( pda::create_problem_eigen(declval<mesh_t>(), declval<prob_t>(), declval<order_t>()) );
      using state_t        = typename app_t::state_type;
      using jacob_t        = typename app_t::jacobian_type;
      using stepper_t      = decltype( pode::create_implicit_stepper(declval<scheme_t>(), declval<app_t &>()) );
      // TODO: generalize
      using linsolver_t    = pls::Solver<pls::iterative::Bicgstab, jacob_t>;
      using nonlinsolver_t = decltype( pnls::create_newton_raphson( declval<stepper_t &>(), declval<linsolver_t&>()) );

    public:

        SchwarzDecomp(
          prob_t probId,
          order_t order,
          scheme_t scheme,
          const std::string & meshRoot,
          vector<double> & dtVec,
          const int icflag = 1)
        {

          // get decomposition info
          read_domain_info(meshRoot);
          ndomains = ndomX * ndomY * ndomZ;

          // set up problem
          setup_controller(dtVec);
          init_problem(probId, order, scheme, meshRoot, icflag);

          // set up communication patterns
          bcStencilSize = (pda::reconstructionTypeToStencilSize(order) - 1) / 2;
          exchDomIdVec = calc_neighbor_dims();
          check_mesh_compat(); // a little error checking
          exchGraphVec = calc_exch_graph(bcStencilSize, exchDomIdVec);
          stateBcVec = init_schw_bc_state();

          // first communication, set internal pointer
          for (int domIdx = 0; domIdx < ndomains; ++domIdx) {
            broadcast_bcState(domIdx);
            appVec[domIdx].setStateBc(&stateBcVec[domIdx]);
          }

          // set up ghost filling graph, set internal pointer
          ghostGraphVec = calc_ghost_graph();
          for (int domIdx = 0; domIdx < ndomains; ++domIdx) {
            // cerr << "Domain " << domIdx << endl;
            // cerr << ghostGraphVec[domIdx] << endl;
            appVec[domIdx].setGraphBc(&ghostGraphVec[domIdx]);
          }

        }

    private:

      void read_domain_info(const string & meshRoot)
      {
        const auto inFile = meshRoot + "/info_domain.dat";
        std::ifstream foundFile(inFile);
        if(!foundFile){
          std::cout << "file not found " << inFile << std::endl;
          exit(EXIT_FAILURE);
        }

        // defaults
        ndomX = 1;
        ndomY = 1;
        ndomZ = 1;

        std::ifstream source( inFile, std::ios_base::in);
        std::string line;
        while (std::getline(source, line) )
        {
          std::istringstream ss(line);
          std::string colVal;
          ss >> colVal;

          if (colVal == "dim"){
            ss >> colVal;
            dim = std::stoi(colVal);
          }

          else if (colVal == "ndomX"){
            ss >> colVal;
            ndomX = std::stoi(colVal);
          }

          else if (colVal == "ndomY"){
            ss >> colVal;
            ndomY = std::stoi(colVal);
          }

          else if (colVal == "ndomZ"){
            ss >> colVal;
            ndomZ = std::stoi(colVal);
          }

          else if (colVal == "overlap"){
            ss >> colVal;
            overlap = std::stoi(colVal);
            // has to be an even number for simplicity, can change later
            if (overlap % 2) {
              std::cerr << "overlap must be an even number" << std::endl;
              exit(-1);
            }
          }
        }
        source.close();
      }

      void setup_controller(vector<double> & dtVec) {

        // physical time step checks
        dt = dtVec;
        if (dt.size() == 1) {
          dt.resize(ndomains, dt[0]);
        } else {
          if (dt.size() != ndomains) {
            cerr << "dt.size() must be 1 or ndomains, exiting" << endl;
            exit(-1);
          }
        }
        dtMax = *max_element(dt.begin(), dt.end());

        // controller time step checks
        controlIters.resize(ndomains);
        for (int domIdx = 0; domIdx < dt.size(); ++domIdx) {
          double niters = dtMax / dt[domIdx];
          if (round(niters) == niters) {
            controlIters[domIdx] = int(round(niters));
          } else {
            cerr << "dt of domain " << domIdx << " (" << dt[domIdx] << ") is not an integer divisor of maximum dt (" << dtMax << ")" << endl;
            exit(-1);
          }
        }

      }

      void init_problem(
        prob_t probId,
        order_t order,
        scheme_t scheme,
        const string & meshRoot,
        const int icflag)
      {
        // problem vectors initialization for each subdomain

        linSolverObj = new linsolver_t;
        meshVec.resize(ndomains);
        appVec.resize(ndomains);
        stateVec.resize(ndomains);
        stateHistVec.resize(ndomains);
        for (int domIdx = 0; domIdx < ndomains; ++domIdx)
        {

          // mesh
          meshVec[domIdx]  = pda::load_cellcentered_uniform_mesh_eigen(meshRoot + "/domain_" + to_string(domIdx));

          // problem and state
          appVec[domIdx]   = pda::create_problem_eigen(meshVec[domIdx], probId, order, icflag);
          stateVec[domIdx] = appVec[domIdx].initialCondition();
          for (int histIdx = 0; histIdx < controlIters[domIdx] + 1; ++histIdx) {
            stateHistVec[domIdx].emplace_back(appVec[domIdx].initialCondition());
          }

          // time stepping
          stepperVec.emplace_back( pode::create_implicit_stepper(scheme, appVec[domIdx]) );
          nonlinSolverVec.emplace_back( pnls::create_newton_raphson(stepperVec.back(), *linSolverObj) );
          nonlinSolverVec[domIdx].setTolerance(1e-5);

        }
      }

      vector<vector<int>> calc_neighbor_dims()
      {
        // determine neighboring domain IDs

        int maxDomNeighbors = 2 * dim;
        vector<vector<int>> exchDomIds(ndomains, vector<int>(maxDomNeighbors, -1));

        for (int domIdx = 0; domIdx < ndomains; ++domIdx) {

          // subdomain indices
          int i, j, k;
          i = domIdx % ndomX;
          if (dim > 1) {
            j = domIdx / ndomX;
          }
          if (dim > 2) {
            k = domIdx / (ndomX * ndomY);
          }


          // 1D, 2D, and 3D
          // left boundary
          if (i != 0) {
            exchDomIds[domIdx][0] = domIdx - 1;
          }

          // right boundary
          if (i != (ndomX - 1)) {
            // ordering change for 1D vs. 2D/3D faces
            if (dim == 1) {
              exchDomIds[domIdx][1] = domIdx + 1;
            }
            else {
              exchDomIds[domIdx][2] = domIdx + 1;
            }
          }

          // 2D and 3D
          if (dim > 1) {
            // front boundary
            if (j != (ndomY - 1)) {
              exchDomIds[domIdx][1] = domIdx + ndomX;
            }

            // back boundary
            if (j != 0) {
              exchDomIds[domIdx][3] = domIdx - ndomX;
            }
          }

          // 3D
          if (dim > 2) {
            // bottom boundary
            if (k != 0) {
              exchDomIds[domIdx][4] = domIdx - (ndomX * ndomY);
            }

            // top boundary
            if (k != (ndomZ - 1)) {
              exchDomIds[domIdx][5] = domIdx + (ndomX * ndomY);
            }
          }

        }

        return exchDomIds;

      }

      void check_mesh_compat() {

        // TODO: extend this for differing (but aligned) mesh resolutions
        if (dim == 1) return; // TODO: still need to check for differing 1D resolutions

        for (int domIdx = 0; domIdx < ndomains; ++domIdx) {

          const auto domMesh = appVec[domIdx].getMesh();
          int nx = domMesh.nx();
          int ny = domMesh.ny();
          int nz = domMesh.nz();

          for (int neighIdx = 0; neighIdx < exchDomIdVec[domIdx].size(); ++neighIdx) {

            int neighDomIdx = exchDomIdVec[domIdx][neighIdx];
            if (neighDomIdx == -1) {
              continue;  // not a Schwarz BC
            }

            const auto neighMesh = appVec[neighDomIdx].getMesh();
            int nxNeigh = neighMesh.nx();
            int nyNeigh = neighMesh.ny();
            int nzNeigh = neighMesh.nz();

            string xerr = "Mesh x-dimension mismatch for domains " + to_string(domIdx) + " v " + to_string(neighDomIdx) + ": " + to_string(nx) + " != " + to_string(nxNeigh);
            string yerr = "Mesh y-dimension mismatch for domains " + to_string(domIdx) + " v " + to_string(neighDomIdx) + ": " + to_string(ny) + " != " + to_string(nyNeigh);
            string zerr = "Mesh z-dimension mismatch for domains " + to_string(domIdx) + " v " + to_string(neighDomIdx) + ": " + to_string(nz) + " != " + to_string(nzNeigh);

            // left and right
            if ((neighIdx == 0) || (neighIdx == 2)) {
              if (ny != nyNeigh) {
                cerr << yerr << endl;
                exit(-1);
              }
              if (nz != nzNeigh) {
                cerr << zerr << endl;
                exit(-1);
              }
            }

            // front and back
            if ((neighIdx == 1) || (neighIdx == 3)) {
              if (nx != nxNeigh) {
                cerr << xerr << endl;
                exit(-1);
              }
              if (nz != nzNeigh) {
                cerr << zerr << endl;
                exit(-1);
              }
            }

            // bottom and top
            if ((neighIdx == 4) || (neighIdx == 5)) {
              if (nx != nxNeigh) {
                cerr << xerr << endl;
                exit(-1);
              }
              if (ny != nyNeigh) {
                cerr << yerr << endl;
                exit(-1);
              }
            }

          } // domain loop
        } // neightbor loop

      }

      vector<graph_t> calc_exch_graph(
        const int bcStencil,
        const vector<vector<int>> & exchDomIds)
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

        vector<graph_t> exchGraphs(ndomains);

        for (int domIdx = 0; domIdx < ndomains; ++domIdx) {

          // this domain's mesh and dimensions
          const auto domMesh = appVec[domIdx].getMesh();
          const auto domGraph = domMesh.graph();
          int nx = domMesh.nx();
          int ny = domMesh.ny();

          const auto & x = domMesh.viewX();
          const auto & y = domMesh.viewY();

          // TODO: generalize to 3D
          pda::resize(exchGraphs[domIdx], 2*nx + 2*ny, bcStencil);
          exchGraphs[domIdx].fill(-1);

          // loop through neighboring domains
          for (int neighIdx = 0; neighIdx < exchDomIds[domIdx].size(); ++neighIdx) {

            int neighDomIdx = exchDomIds[domIdx][neighIdx];
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
              int bcCellIdx = 0; // left boundary is the start
              for (int yIdx = 0; yIdx < ny; ++yIdx) {
                for (int stencilIdx = 0; stencilIdx < bcStencil; ++stencilIdx) {
                  // exchCellIdx = (nxNeigh * (yIdx + 1) - 1) - overlap - stencilIdx;
                  exchCellIdx = (nxNeigh * (yIdx + 1) - 1) - overlap + 1 - stencilIdx; // toward boundary
                  exchGraphs[domIdx](bcCellIdx, stencilIdx) = exchCellIdx;
                }
                bcCellIdx++;
              }
            }

            // north-south neighbors will share a column index
            // "front"
            if (neighIdx == 1) {
              int bcCellIdx = ny * bcStencil;  // skip left boundary indices
              for (int xIdx = 0; xIdx < nx; ++xIdx) {
                for (int stencilIdx = 0; stencilIdx < bcStencil; ++stencilIdx) {
                  exchCellIdx = (overlap + stencilIdx) * nxNeigh + xIdx;
                  exchGraphs[domIdx](bcCellIdx, stencilIdx) = exchCellIdx;
                }
                bcCellIdx++;
              }
            }

            // right
            if (neighIdx == 2) {
              int bcCellIdx = (ny + nx) * bcStencil; // skip left and "front" boundary indices
              for (int yIdx = 0; yIdx < ny; ++yIdx) {
                for (int stencilIdx = 0; stencilIdx < bcStencil; ++stencilIdx) {
                  exchCellIdx = nxNeigh * yIdx + overlap + stencilIdx;
                  // TODO: check for a different problem
                  // exchCellIdx = nxNeigh * yIdx + overlap - 1 + stencilIdx; // toward boundary
                  // exchCellIdx = nxNeigh * yIdx + overlap + 1 + stencilIdx; // away from boundary
                  exchGraphs[domIdx](bcCellIdx, stencilIdx) = exchCellIdx;
                }
                bcCellIdx++;
              }
            }

            // "back"
            if (neighIdx == 3) {
              int bcCellIdx = (2*ny + nx) * bcStencil;  // skip left, "front", and right boundary indices
              for (int xIdx = 0; xIdx < nx; ++xIdx) {
                for (int stencilIdx = 0; stencilIdx < bcStencil; ++stencilIdx) {
                  exchCellIdx = (nyNeigh - 1 - overlap - stencilIdx) * nxNeigh + xIdx;
                  exchGraphs[domIdx](bcCellIdx, stencilIdx) = exchCellIdx;
                }
                bcCellIdx++;
              }
            }

            // TODO: generalize to 3D

          } // neighbor loop
        } // domain loop

        return exchGraphs;

      }

      vector<state_t> init_schw_bc_state() {

        int bcStencilDof = bcStencilSize * app_t::numDofPerCell;
        vector<state_t> stateBcs(ndomains);

        for (int domIdx = 0; domIdx < ndomains; ++domIdx) {
          const auto meshObj = meshVec[domIdx];
          int numDofStencilBc = 2 * bcStencilDof * (meshObj.nx() + meshObj.ny() + meshObj.nz());
          pda::resize(stateBcs[domIdx], numDofStencilBc);
          stateBcs[domIdx].fill(0.0);
        }

        return stateBcs;

      }

      void comm_stateBc(
        const int startIdx,
        const int endIdx,
        const graph_t & exchGraph,
        state_t & bcState,
        const state_t & intState)
      {

        int dofPerCell = app_t::numDofPerCell;
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

      void broadcast_bcState(const int domIdx)
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
          int nzNeigh = neighMesh.nz();
          auto * neighStateBc = &stateBcVec[neighDomIdx];
          const auto neighExchGraph = exchGraphVec[neighDomIdx];

          // TODO: extend to 3D, need to change L/R and F/B indices to account for nzNeigh

          // this domain is the neighboring domain's left neighbor
          if (neighIdx == 2) {
            startIdx = 0;
            endIdx = nyNeigh;
            comm_stateBc(startIdx, endIdx, neighExchGraph, *neighStateBc, domState);
          }

          // this domain is the neighboring domain's front neighbor
          if (neighIdx == 3) {
            startIdx = nyNeigh;
            endIdx = nyNeigh + nxNeigh;
            comm_stateBc(startIdx, endIdx, neighExchGraph, *neighStateBc, domState);
          }

          // this domain is the neighboring domain's right neighbor
          if (neighIdx == 0) {
            startIdx = nyNeigh + nxNeigh;
            endIdx = 2 * nyNeigh + nxNeigh;
            comm_stateBc(startIdx, endIdx, neighExchGraph, *neighStateBc, domState);
          }

          // this domain is the neighboring domain's back neighbor
          if (neighIdx == 1) {
            startIdx = 2 * nyNeigh + nxNeigh;
            endIdx = 2 * nyNeigh + 2 * nxNeigh;
            comm_stateBc(startIdx, endIdx, neighExchGraph, *neighStateBc, domState);
          }

          // this domain is the neighboring domain's bottom neighbor
          // if (neighIdx == 5) {
          //   startIdx = ;
          //   endIdx = ;
          // }

          // this domain is the neighboring domain's top neighbor
          // if (neighIdx == 4) {
          //   startIdx = ;
          //   endIdx = ;
          // }

          // comm_stateBc(startIdx, endIdx, neighExchGraph, *neighStateBc, domState);

        }

      }

      vector<graph_t> calc_ghost_graph()
      {

        int dofPerCell = app_t::numDofPerCell;
        vector<graph_t> ghostGraphs(ndomains);

        for (int domIdx = 0; domIdx < ndomains; ++domIdx) {

          const auto meshObj = meshVec[domIdx];
          const auto intGraph = meshObj.graph();
          int nx = meshObj.nx();
          int ny = meshObj.ny();
          int nz = meshObj.nz();

          const auto & rowsBd = meshObj.graphRowsOfCellsNearBd();
          pda::resize(ghostGraphs[domIdx], int(rowsBd.size()), 2 * dim);
          ghostGraphs[domIdx].fill(-1);

          for (decltype(rowsBd.size()) it = 0; it < rowsBd.size(); ++it) {

            // ASK FR: is there an instance when rowsBd[it] != intGraph(rowsBd[it], 0)?
            //    The indices appear to be identical
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
                ghostGraphs[domIdx](it, 0) = bcCellIdx * dofPerCell;
              }
            }

            if (front0 == -1) {
              if (exchDomIdVec[domIdx][1] != -1) {
                bcCellIdx = ny + colIdx;
                ghostGraphs[domIdx](it, 1) = bcCellIdx * dofPerCell;
              }
            }

            if (right0 == -1) {
              if (exchDomIdVec[domIdx][2] != -1) {
                bcCellIdx = ny + nx + rowIdx;
                ghostGraphs[domIdx](it, 2) = bcCellIdx * dofPerCell;
              }
            }

            if (back0 == -1) {
              if (exchDomIdVec[domIdx][3] != -1) {
                bcCellIdx = 2 * ny + nx + colIdx;
                ghostGraphs[domIdx](it, 3) = bcCellIdx * dofPerCell;
              }
            }
            // TODO: extend to higher order, 3D

          } // boundary cell loop
        } // domain loop

        return ghostGraphs;

      }

      array<double, 2> calcConvergence(const state_t & state1, const state_t & state2)
      {
        // TODO: assumed to be an Eigen state, not sure how to generalize
        // TODO: compute convergence for each variable separately

        int numDOF = state1.size();
        if (state2.size() != numDOF) {
          cerr << "state1 size does not match state2 size, " << numDOF << " vs. " << state2.size() << endl;
          exit(-1);
        }

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

    public:

      void calc_controller_step(
        int outerStep,
        double time,
        const double rel_err_tol,
        const double abs_err_tol,
        const int convergeStepMax)
      {

        // store initial step for resetting if Schwarz iter does not converge
        for (int domIdx = 0; domIdx < ndomains; ++domIdx) {
          stateHistVec[domIdx][0] = stateVec[domIdx];
        }

        // convergence
        int convergeStep = 0;
        vector<array<double, 2>> convergeVals(ndomains);
        while (convergeStep < convergeStepMax) {

          cout << "Schwarz iteration " << convergeStep + 1 << endl;

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
                convergeVals[domIdx] = calcConvergence(stateVec[domIdx], stateHistVec[domIdx].back());
              }

              // store intra-step history
              stateHistVec[domIdx][innerStep + 1] = stateVec[domIdx];

              // set (interpolated) boundary conditions

              // update local step and time
              stepDom++;
              timeDom += dtDom;

            } // domain loop

            // broadcast boundary conditions
            broadcast_bcState(domIdx);

          }

          // check convergence for all domains, break if conditions met
          double abs_err = 0.0;
          double rel_err = 0.0;
          for (int domIdx = 0; domIdx < ndomains; ++domIdx) {
            abs_err += convergeVals[domIdx][0];
            rel_err += convergeVals[domIdx][1];
          }
          abs_err /= ndomains;
          rel_err /= ndomains;
          cout << "Average abs err: " << abs_err << endl;
          cout << "Average rel err: " << rel_err << endl;
          if ((rel_err < rel_err_tol) || (abs_err < abs_err_tol)) {
            break;
          }

          convergeStep++;

          // reset interior state if not converged
          for (int domIdx = 0; domIdx < ndomains; ++domIdx) {
            stateVec[domIdx] = stateHistVec[domIdx][0];
          }

        } // convergence loop

      }


    public:

      int ndomains;
      double dtMax;
      vector<state_t> stateVec;

    private:

      // mesh decomposition
      int dim;
      int ndomX = 1;
      int ndomY = 1;
      int ndomZ = 1;
      int overlap;

      // subdomain communication
      int bcStencilSize;
      int bcStencilDof;
      vector<vector<int>> exchDomIdVec;
      vector<graph_t> exchGraphVec;
      vector<state_t> stateBcVec;
      vector<graph_t> ghostGraphVec;

      // time-stepping
      vector<double> dt;
      vector<int> controlIters;

      // problem vectors
      linsolver_t* linSolverObj;
      vector<mesh_t> meshVec;
      vector<app_t> appVec;
      vector<vector<state_t>> stateHistVec;
      vector<stepper_t> stepperVec;
      vector<nonlinsolver_t> nonlinSolverVec;

  };
}}

#endif