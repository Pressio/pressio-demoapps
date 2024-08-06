/*
//@HEADER
// ************************************************************************
//
// advection_diffusion_2d_prob_class.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef PRESSIODEMOAPPS_ADVECTION_DIFFUSION_2D_PROB_CLASS_HPP_
#define PRESSIODEMOAPPS_ADVECTION_DIFFUSION_2D_PROB_CLASS_HPP_

#include "advection_diffusion_2d_flux_functions.hpp"
#include "advection_diffusion_2d_initial_condition.hpp"
#include "advection_diffusion_2d_ghost_filler_neumann.hpp"
#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "advection_diffusion_2d_flux_mixin.hpp"
#include "mixin_directional_flux_balance.hpp"
#include "mixin_directional_flux_balance_jacobian.hpp"
#include "Eigen/Sparse"
#include "noop.hpp"
#include "custom_bcs_functions.hpp"
#include "ghost_relative_locations.hpp"

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#include <omp.h>
#endif

namespace pressiodemoapps{ namespace impladvdiff2d{

// tags are used inside he public create function: create_problem_...()
// in the file ../advection_diffusion.hpp
// to dispatch to the proper problem
// so add new ones if a new problem is added
struct TagBurgersPeriodic{};


/////////////////////////
// eigen class
/////////////////////////
template<
  class MeshType,
  class BCFunctorsHolderType = impl::NoOperation<void>
  >
class EigenApp
{

public:
  // required public aliases
  using index_t	      = typename MeshType::index_t;
  using scalar_type   = typename MeshType::scalar_t;
  using state_type    = Eigen::Matrix<scalar_type,Eigen::Dynamic,1>;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, index_t>;
  using mesh_connectivity_graph_type = typename MeshType::graph_t;

  static constexpr int dimensionality{2};
  static constexpr int numDofPerCell{2};

private:

  using ghost_container_type      = Eigen::Matrix<scalar_type,
						  Eigen::Dynamic,
						  Eigen::Dynamic,
						  Eigen::RowMajor>;
  using stencil_container_type    = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using flux_type	          = Eigen::Matrix<scalar_type, numDofPerCell, 1>;
  using edge_rec_type	          = Eigen::Matrix<scalar_type, numDofPerCell, 1>;
  using flux_jac_type             = Eigen::Matrix<scalar_type, numDofPerCell, numDofPerCell>;
  using reconstruction_gradient_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  EigenApp() = delete;

  //
  // constructor for Burgers2d
  //
  EigenApp(const MeshType & meshObj,
	   ::pressiodemoapps::AdvectionDiffusion2d probEnum,
	   ::pressiodemoapps::InviscidFluxReconstruction inviscidFluxRecEn,
	   ::pressiodemoapps::InviscidFluxScheme invFluxSchemeEn,
	   ::pressiodemoapps::ViscousFluxReconstruction visFluxRecEn,
	   scalar_type icPulseMagnitude,
	   scalar_type icSpread,
	   scalar_type diffusionCoeff,
	   scalar_type x0,
	   scalar_type y0)
    : m_probEn(probEnum),
      m_inviscidFluxRecEn(inviscidFluxRecEn),
      m_inviscidFluxSchemeEn(invFluxSchemeEn),
      m_viscousFluxRecEn(visFluxRecEn),
      m_meshObj(meshObj),
      m_burgers2d_icPulse(icPulseMagnitude),
      m_burgers2d_icSpread(icSpread),
      m_burgers2d_diffusion(diffusionCoeff),
      m_burgers2d_x0(x0),
      m_burgers2d_y0(y0)
  {
    m_numDofStencilMesh = m_meshObj.get().stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.get().sampleMeshSize()  * numDofPerCell;

    if (m_probEn == pressiodemoapps::AdvectionDiffusion2d::BurgersOutflow) {
      allocateGhosts();
    }

  }

// #if !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
//   EigenApp(const MeshType & meshObj,
// 	   ::pressiodemoapps::AdvectionDiffusion2d probEnum,
// 	   ::pressiodemoapps::InviscidFluxReconstruction inviscidFluxRecEn,
// 	   ::pressiodemoapps::InviscidFluxScheme invFluxSchemeEn,
// 	   ::pressiodemoapps::ViscousFluxReconstruction visFluxRecEn,
// 	   scalar_type icPulseMagnitude,
// 	   scalar_type icSpread,
// 	   scalar_type diffusionCoeff,
// 	   scalar_type x0,
// 	   scalar_type y0,
// 	   BCFunctorsHolderType && bcHolder)
//     : m_probEn(probEnum),
//       m_inviscidFluxRecEn(inviscidFluxRecEn),
//       m_inviscidFluxSchemeEn(invFluxSchemeEn),
//       m_viscousFluxRecEn(visFluxRecEn),
//       m_meshObj(meshObj),
//       m_burgers2d_icPulse(icPulseMagnitude),
//       m_burgers2d_icSpread(icSpread),
//       m_burgers2d_diffusion(diffusionCoeff),
//       m_burgers2d_x0(x0),
//       m_burgers2d_y0(y0),
//       m_bcFuncsHolder(std::move(bcHolder))
//   {
//     m_numDofStencilMesh = m_meshObj.get().stencilMeshSize() * numDofPerCell;
//     m_numDofSampleMesh  = m_meshObj.get().sampleMeshSize()  * numDofPerCell;

//     if (m_probEn == pressiodemoapps::AdvectionDiffusion2d::BurgersOutflow) {
// 	allocateGhosts();
//     }

//   }
// #endif

  state_type initialCondition() const
  {
    state_type initialState(m_numDofStencilMesh);

    if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersPeriodic ||
	m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersOutflow)
      {
	burgers2d_gaussian(initialState, m_meshObj.get(),
			   m_burgers2d_icPulse,
			   m_burgers2d_icSpread,
			   m_burgers2d_x0,
			   m_burgers2d_y0);
      }

    return initialState;
  }

// public:
//   template <class T>
//   void setBCPointer(::pressiodemoapps::impl::GhostRelativeLocation rloc, T* ptr) {
//     m_bcFuncsHolder.setInternalPointer(rloc, ptr);
//   }

protected:
  int numDofPerCellImpl() const {
    return numDofPerCell;
  }

  void initializeJacobian(jacobian_type & J)
  {
    J.resize(m_numDofSampleMesh, m_numDofStencilMesh);

    using Tr = Eigen::Triplet<scalar_type>;
    std::vector<Tr> trList;
    initializeJacobianForNearBoundaryCells(trList);
    initializeJacobianForInnerCells(trList);
    J.setFromTriplets(trList.begin(), trList.end());

    // compress to make it Csr
    if (!J.isCompressed()){
      J.makeCompressed();
    }
  }

  // note that we MUST template U_t, V_t because
  // when doing bindings, these are deduced to be a Eigen Ref
  template<class U_t, class V_t>
  void velocityAndOptionalJacobian(const U_t & U,
				   const scalar_type currentTime,
				   V_t & V,
				   jacobian_type * J) const
  {

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp parallel
    {
#endif

      // for omp, these are private variables for each thread
      // edge reconstructions
      edge_rec_type uMinusHalfNeg, uMinusHalfPos;
      edge_rec_type uPlusHalfNeg, uPlusHalfPos;
      // fluxes
      flux_type fluxL, fluxF;
      flux_type fluxR, fluxB;

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
      ::pressiodemoapps::set_zero_omp(V);
#else
      ::pressiodemoapps::set_zero(V);
#endif

      if (J){
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
	::pressiodemoapps::set_zero_omp(*J);
#else
	::pressiodemoapps::set_zero(*J);
#endif
      }

      if constexpr(std::is_same_v<BCFunctorsHolderType, impl::NoOperation<void>>){
        fillGhosts(U, currentTime);
      }
      else {
	throw std::runtime_error("Custom BCs not implemented yet")
      	// fillGhostsUseCustomFunctors(U, currentTime, m_meshObj, m_bcFuncsHolder,
	// 			    m_ghostLeft, m_ghostFront,
	// 			    m_ghostRight, m_ghostBack, numDofPerCell);
      }

      if (J){
	velocityAndJacobianImpl(U, currentTime, V, *J,
				fluxL, fluxF, fluxR, fluxB,
				uMinusHalfNeg, uMinusHalfPos,
				uPlusHalfNeg,  uPlusHalfPos);
      }
      else{
        if (!m_meshObj.get().isFullyPeriodic()) {
          velocityOnlyNearBdCellsImpl(U, currentTime, V,
				      fluxL, fluxF, fluxR, fluxB,
				      uMinusHalfNeg, uMinusHalfPos,
				      uPlusHalfNeg,  uPlusHalfPos);
	}

	velocityOnlyInnerCellsImpl(U, currentTime, V,
				   fluxL, fluxF, fluxR, fluxB,
				   uMinusHalfNeg, uMinusHalfPos,
				   uPlusHalfNeg,  uPlusHalfPos);
      }

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
}//end omp parallel
#endif

  }

private:
  template<class U_t>
  void fillGhosts(const U_t & U, const scalar_type currentTime) const
  {
    const auto stencilSize = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    if (m_probEn == pressiodemoapps::AdvectionDiffusion2d::BurgersOutflow)
    {
      using ghost_filler_t  = Ghost2dNeumannFiller<
	U_t, MeshType, ghost_container_type>;
      ghost_filler_t ghF(stencilSize, numDofPerCell,
			 U, m_meshObj.get(),
			 m_ghostLeft, m_ghostFront,
			 m_ghostRight, m_ghostBack);
      const auto & rowsBd = m_meshObj.get().graphRowsOfCellsNearBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
      for (decltype(rowsBd.size()) it=0; it<rowsBd.size(); ++it){
	      ghF(rowsBd[it], it);
      }
    }
    else {
      // no op
    }
  }

  template<class Tr>
  void initializeJacobianForInnerCells(std::vector<Tr> & trList)
  {
    // for inner cells, the Jacobian is exact and depends on
    // the scheme wanted by the user, no special treatment needed

    const scalar_type zero = 0;
    const auto & graph = m_meshObj.get().graph();
    // only grab the graph rows for INNER cells (i.e. AWAY from boundaries)
    const auto & targetGraphRows = m_meshObj.get().graphRowsOfCellsAwayFromBd();
    for (std::size_t it=0; it<targetGraphRows.size(); ++it)
      {
	const auto smPt = targetGraphRows[it];

	// find out which row in the jacobian we are dealing with
	const auto jacRowOfCurrCellFirstDof = smPt*numDofPerCell;

	// initialize jacobian block entries wrt current cell's dofs
	const auto jacColOfCurrCellRho = graph(smPt, 0)*numDofPerCell;
	for (int k=0; k<numDofPerCell; ++k){
	  for (int j=0; j<numDofPerCell; ++j){
	    trList.push_back( Tr(jacRowOfCurrCellFirstDof+k,
				 jacColOfCurrCellRho+j,
				 zero) );
	  }
	}

	// now fill wrt neighbors: this depends on the scheme
	// so I need to use the largest between the stencil needed
	// for inviscid and viscous reconstructions.
	// Most likely the inviscid always wins so use that for now.
	// if assert breaks, this will need to be fixed.
	assert(reconstructionTypeToStencilSize(m_inviscidFluxRecEn)
	       >= reconstructionTypeToStencilSize(m_viscousFluxRecEn));

	// find out how many neighboring cells: recall this is 2d
	const int numNeighbors =
	  (m_inviscidFluxRecEn == InviscidFluxReconstruction::FirstOrder) ? 4
	    : (m_inviscidFluxRecEn == InviscidFluxReconstruction::Weno3) ? 8
	      : (m_inviscidFluxRecEn == InviscidFluxReconstruction::Weno5) ? 12 : -1;
	assert(numNeighbors > 0);

	for (int i=1; i<=numNeighbors; ++i){
	  const auto colInd = graph(smPt, i)*numDofPerCell;
	  for (int k=0; k<numDofPerCell; ++k){
	    for (int j=0; j<numDofPerCell; ++j){
	      trList.push_back( Tr(jacRowOfCurrCellFirstDof+k, colInd+j, zero) );
	    }
	  }
	}
      }
  }

  template<class Tr>
  void initializeJacobianForNearBoundaryCells(std::vector<Tr> & trList)
  {
    const scalar_type zero = 0;
    const auto & graph = m_meshObj.get().graph();
    const auto & targetGraphRows = m_meshObj.get().graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<targetGraphRows.size(); ++it)
      {
	const auto smPt = targetGraphRows[it];
	const auto jacRowOfCurrCellRho = smPt*numDofPerCell;
	const auto jacColOfCurrCellRho = graph(smPt, 0)*numDofPerCell;

	// wrt current cell's dofs
	for (int k=0; k<numDofPerCell; ++k){
	  for (int j=0; j<numDofPerCell; ++j){
	    trList.push_back( Tr(jacRowOfCurrCellRho+k, jacColOfCurrCellRho+j, zero) );
	  }
	}

	// wrt neighbors
	// for near-bd, we only do first-order Jacobian for now
	for (int i=1; i<=4; ++i){
	  const auto nID = graph(smPt, i);
	  if (nID != -1){
	    const auto ci = nID*numDofPerCell;
	    for (int k=0; k<numDofPerCell; ++k){
	      for (int j=0; j<numDofPerCell; ++j){
		trList.push_back( Tr(jacRowOfCurrCellRho+k, ci+j, zero) );
	      }
	    }
	  }
	}
      }
  }


  template<class U_t, class V_t>
  void velocityAndJacobianImpl(const U_t & U,
				       const scalar_type currentTime,
				       V_t & V,
				       jacobian_type & J,
				       flux_type & fluxL,
				       flux_type & fluxF,
				       flux_type & fluxR,
				       flux_type & fluxB,
				       edge_rec_type & uMinusHalfNeg,
				       edge_rec_type & uMinusHalfPos,
				       edge_rec_type & uPlusHalfNeg,
				       edge_rec_type & uPlusHalfPos) const
  {

    // flux jacobians
    flux_jac_type fluxJacLNeg;
    flux_jac_type fluxJacLPos;
    flux_jac_type fluxJacFNeg;
    flux_jac_type fluxJacFPos;
    flux_jac_type fluxJacRNeg;
    flux_jac_type fluxJacRPos;
    flux_jac_type fluxJacBNeg;
    flux_jac_type fluxJacBPos;

    int nonZerosCountBeforeComputing = J.nonZeros();

    if (!m_meshObj.get().isFullyPeriodic()){
      if (m_inviscidFluxRecEn == InviscidFluxReconstruction::FirstOrder){
	velocityAndJacNearBDCellsImplFirstOrder(U, currentTime, V, J,
						fluxL, fluxF, fluxR, fluxB,
						fluxJacLNeg, fluxJacLPos,
						fluxJacFNeg, fluxJacFPos,
						fluxJacRNeg, fluxJacRPos,
						fluxJacBNeg, fluxJacBPos,
						uMinusHalfNeg, uMinusHalfPos,
						uPlusHalfNeg,  uPlusHalfPos);
      }
      else{
	velocityAndJacNearBDCellsImplDifferentSchemeNotFirstOrderInviscid(U, currentTime, V, J,
									  fluxL, fluxF, fluxR, fluxB,
									  fluxJacLNeg, fluxJacLPos,
									  fluxJacFNeg, fluxJacFPos,
									  fluxJacRNeg, fluxJacRPos,
									  fluxJacBNeg, fluxJacBPos,
									  uMinusHalfNeg, uMinusHalfPos,
									  uPlusHalfNeg,  uPlusHalfPos);
      }
    }

    velocityAndJacInnerCellsImpl(U, currentTime, V, J,
				 fluxL, fluxF, fluxR, fluxB,
				 fluxJacLNeg, fluxJacLPos,
				 fluxJacFNeg, fluxJacFPos,
				 fluxJacRNeg, fluxJacRPos,
				 fluxJacBNeg, fluxJacBPos,
				 uMinusHalfNeg, uMinusHalfPos,
				 uPlusHalfNeg,  uPlusHalfPos);

    assert(J.nonZeros() == nonZerosCountBeforeComputing);
    (void) nonZerosCountBeforeComputing;
  }

  template<class U_t, class V_t>
  void velocityAndJacInnerCellsImpl(const U_t & U,
				    const scalar_type /*currentTime*/,
				    V_t & V,
				    jacobian_type & J,
				    flux_type & fluxL,
				    flux_type & fluxF,
				    flux_type & fluxR,
				    flux_type & fluxB,
				    flux_jac_type & fluxJacLNeg,
				    flux_jac_type & fluxJacLPos,
				    flux_jac_type & fluxJacFNeg,
				    flux_jac_type & fluxJacFPos,
				    flux_jac_type & fluxJacRNeg,
				    flux_jac_type & fluxJacRPos,
				    flux_jac_type & fluxJacBNeg,
				    flux_jac_type & fluxJacBPos,
				    edge_rec_type & uMinusHalfNeg,
				    edge_rec_type & uMinusHalfPos,
				    edge_rec_type & uPlusHalfNeg,
				    edge_rec_type & uPlusHalfPos) const
  {
    // for inner cells, velocity and Jacobian
    // are both computed according to the order selected by the user
    // because for inner cells we support also Jacobians for Weno

    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;

    // allocate gradients of reconstructed states
    // the size depends on the scheme selected
    const auto stencilSize = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    assert(stencilSize >= reconstructionTypeToStencilSize(m_viscousFluxRecEn));

    reconstruction_gradient_t gradLNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradLPos(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradFNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradFPos(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRPos(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradBNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradBPos(numDofPerCell, stencilSize-1);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impl::ComputeDirectionalFluxBalanceJacobianOnInteriorCell<
	  pda::impladvdiff2d::ComputeDirectionalFluxValuesAndJacobians<
	    pda::impl::ReconstructorForDiscreteFunction<
	      dimensionality, MeshType, U_t, edge_rec_type, reconstruction_gradient_t>,
	    scalar_type, flux_type, flux_jac_type>,
	  dimensionality, MeshType, jacobian_type>,
      V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.get().dxInv(),
		    /* end args for velo */
		    J, xAxis, m_meshObj.get(),
		    /* end args for jac */
		    m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR,
		    fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
		    /* end args for flux */
		    xAxis, toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj.get(),
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		    gradLNeg, gradLPos, gradRNeg, gradRPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.get().dyInv(),
		    /* end args for velo */
		    J, yAxis, m_meshObj.get(),
		    /* end args for jac */
		    m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF,
		    fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
		    /* end args for flux */
		    yAxis, toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj.get(),
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		    gradBNeg, gradBPos, gradFNeg, gradFPos
		    /* end args for reconstructor */
		    );

    const auto & graphRows = m_meshObj.get().graphRowsOfCellsAwayFromBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];
      Fx(smPt, numDofPerCell);
      Fy(smPt, numDofPerCell);

      // diffusion contribution
      addBurgersDiffusionToVelocityAndOptionalJacobianInnerCells(U, V, &J, smPt);
    }
  }

  template<class state_t, class V_t>
  void velocityAndJacNearBDCellsImplDifferentSchemeNotFirstOrderInviscid(
						    const state_t & state,
						    const scalar_type currentTime,
						    V_t & V,
						    jacobian_type & J,
						    flux_type & fluxL,
						    flux_type & fluxF,
						    flux_type & fluxR,
						    flux_type & fluxB,
						    flux_jac_type & fluxJacLNeg,
						    flux_jac_type & fluxJacLPos,
						    flux_jac_type & fluxJacFNeg,
						    flux_jac_type & fluxJacFPos,
						    flux_jac_type & fluxJacRNeg,
						    flux_jac_type & fluxJacRPos,
						    flux_jac_type & fluxJacBNeg,
						    flux_jac_type & fluxJacBPos,
						    edge_rec_type & uMinusHalfNeg,
						    edge_rec_type & uMinusHalfPos,
						    edge_rec_type & uPlusHalfNeg,
						    edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;

    throw std::runtime_error("Higher order not fixed");

    // if here, then the velocity must be computed with Weno,
    /// while the jacobian must be computed with first order

    using stencil_filler_t  = pda::impl::StencilFiller<dimensionality, stencil_container_type,
						       state_t, MeshType, ghost_container_type>;

    // *****************************
    // *** functors for velocity ***
    // *****************************
    const auto stencilSizeForV = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    // this method should only be called for stencilsize 5 or 7, because
    // the case =3 is handled specifically in another method
    assert(stencilSizeForV == 5 or stencilSizeForV == 7);

    stencil_container_type stencilValsForV(numDofPerCell*stencilSizeForV);

    stencil_filler_t FillStencilVeloX(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				      state, m_meshObj.get(), m_ghostLeft, m_ghostRight,
				      stencilValsForV, xAxis);
    stencil_filler_t FillStencilVeloY(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				      state, m_meshObj.get(), m_ghostBack, m_ghostFront,
				      stencilValsForV, yAxis);

    using velo_functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impladvdiff2d::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  scalar_type, flux_type>,
      V_t, scalar_type
      >;

    velo_functor_type funcVeloX(V, m_meshObj.get().dxInv(),
				/* end args for velo */
				m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR,
				/* end args for flux */
				toReconstructionScheme(m_inviscidFluxRecEn), stencilValsForV,
				uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
				/* end args for reconstructor */
				);

    velo_functor_type funcVeloY(V, m_meshObj.get().dyInv(),
				/* end args for velo */
				m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF,
				/* end args for flux */
				toReconstructionScheme(m_inviscidFluxRecEn), stencilValsForV,
				uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
				/* end args for reconstructor */
				);

    // *****************************
    // *** functors for jacobian ***
    // *****************************
    const auto firstOrderRec = pda::InviscidFluxReconstruction::FirstOrder;
    const auto stencilSizeForJ = reconstructionTypeToStencilSize(firstOrderRec);
    stencil_container_type stencilValsForJ(numDofPerCell*stencilSizeForJ);
    stencil_filler_t FillStencilJacX(stencilSizeForJ,
				     state, m_meshObj.get(), m_ghostLeft, m_ghostRight,
				     stencilValsForJ, xAxis);
    stencil_filler_t FillStencilJacY(stencilSizeForJ,
				     state, m_meshObj.get(), m_ghostBack, m_ghostFront,
				     stencilValsForJ, yAxis);

    using jac_functor_type =
      pda::impl::ComputeDirectionalFluxBalanceFirstOrderJacobianOnBoundaryCell<
	pda::impladvdiff2d::ComputeDirectionalFluxJacobians<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  scalar_type, flux_jac_type>,
      dimensionality, MeshType, jacobian_type
      >;

    jac_functor_type funcJacX(J, xAxis, m_meshObj.get(),
			      /* end args for jac */
			      m_probEn, m_inviscidFluxSchemeEn, normalX_,
			      fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
			      /* end args for flux */
			      toReconstructionScheme(firstOrderRec), stencilValsForJ,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			      /* end args for reconstructor */
			      );

    jac_functor_type funcJacY(J, yAxis, m_meshObj.get(),
			      /* end args for jac */
			      m_probEn, m_inviscidFluxSchemeEn, normalY_,
			      fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
			      /* end args for flux */
			      toReconstructionScheme(firstOrderRec), stencilValsForJ,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			      /* end args for reconstructor */
			      );

    // ************
    // loop
    // ************
    const auto & graph     = m_meshObj.get().graph();
    const auto & graphRows = m_meshObj.get().graphRowsOfCellsNearBd();
    const auto dxInvSq	   = m_meshObj.get().dxInv()*m_meshObj.get().dxInv();
    const auto dyInvSq	   = m_meshObj.get().dyInv()*m_meshObj.get().dyInv();
    const auto diffDxInvSq = m_burgers2d_diffusion*dxInvSq;
    const auto diffDyInvSq = m_burgers2d_diffusion*dyInvSq;
    constexpr auto two      = static_cast<scalar_type>(2);

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];
      const auto uIndex	= graph(smPt, 0);
      const auto uIndexLeft  = graph(smPt, 1);
      const auto uIndexFront = graph(smPt, 2);
      const auto uIndexRight = graph(smPt, 3);
      const auto uIndexBack  = graph(smPt, 4);

      // x-direction
      FillStencilVeloX(smPt, it, numDofPerCell);
      funcVeloX(smPt, numDofPerCell);
      FillStencilJacX(smPt, it, numDofPerCell);
      if constexpr(std::is_same_v<BCFunctorsHolderType, impl::NoOperation<void>>){
        fillJacFactorsForCellBd(smPt, xAxis);
      }
      else{
        throw std::runtime_error("Need to fix higher order custom BCs");
      }
      funcJacX(smPt, numDofPerCell, m_bcCellJacFactors);
      // diffusion contribution to velocity
      if (stencilSizeForV == 5){
	V(smPt) += diffDxInvSq*( stencilValsForV(3) - two*stencilValsForV(2) + stencilValsForV(1) );
      }
      else if (stencilSizeForV == 7){
	V(smPt) += diffDxInvSq*( stencilValsForV(4) - two*stencilValsForV(3) + stencilValsForV(2) );
      }

      // y-direction
      FillStencilVeloY(smPt, it, numDofPerCell);
      funcVeloY(smPt, numDofPerCell);
      FillStencilJacY(smPt, it, numDofPerCell);
      if constexpr(std::is_same_v<BCFunctorsHolderType, impl::NoOperation<void>>){
        fillJacFactorsForCellBd(smPt, yAxis);
      }
      else{
        throw std::runtime_error("Need to fix higher order custom BCs");
      }
      funcJacY(smPt, numDofPerCell, m_bcCellJacFactors);
      // diffusion contribution to velocity
      if (stencilSizeForV == 5){
	V(smPt) += diffDyInvSq*( stencilValsForV(3) -two*stencilValsForV(2) +stencilValsForV(1) );
      }
      else if (stencilSizeForV == 7){
	V(smPt) += diffDyInvSq*( stencilValsForV(4) -two*stencilValsForV(3) +stencilValsForV(2) );
      }

      // diffusion contribution to Jacobian
      auto selfValue = -two*diffDxInvSq -two*diffDyInvSq;
      if (uIndexLeft != -1){
	J.coeffRef(smPt, uIndexLeft) += diffDxInvSq;
      }
      else{
	selfValue += -diffDxInvSq;
      }

      if (uIndexFront != -1){
	J.coeffRef(smPt, uIndexFront) += diffDyInvSq;
      }else{
	selfValue += -diffDyInvSq;
      }

      if (uIndexRight != -1){
	J.coeffRef(smPt, uIndexRight) += diffDxInvSq;
      }
      else{
	selfValue += -diffDxInvSq;
      }

      if (uIndexBack != -1){
	J.coeffRef(smPt, uIndexBack) += diffDyInvSq;
      }
      else{
	selfValue += -diffDyInvSq;
      }
      J.coeffRef(smPt, uIndex) += selfValue;
    }
  }

  template<class state_t, class V_t>
  void velocityAndJacNearBDCellsImplFirstOrder(const state_t & state,
					       const scalar_type currentTime,
					       V_t & V,
					       jacobian_type & J,
					       flux_type & fluxL,
					       flux_type & fluxF,
					       flux_type & fluxR,
					       flux_type & fluxB,
					       flux_jac_type & fluxJacLNeg,
					       flux_jac_type & fluxJacLPos,
					       flux_jac_type & fluxJacFNeg,
					       flux_jac_type & fluxJacFPos,
					       flux_jac_type & fluxJacRNeg,
					       flux_jac_type & fluxJacRPos,
					       flux_jac_type & fluxJacBNeg,
					       flux_jac_type & fluxJacBPos,
					       edge_rec_type & uMinusHalfNeg,
					       edge_rec_type & uMinusHalfPos,
					       edge_rec_type & uPlusHalfNeg,
					       edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    assert(m_inviscidFluxRecEn == InviscidFluxReconstruction::FirstOrder);

    // if here, then the scheme for velocity matches
    // the one for Jacobian so we can use same functors

    const auto stencilSize = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    stencil_container_type stencilVals(numDofPerCell*stencilSize);

    using stencil_filler_t  = pda::impl::StencilFiller<
    dimensionality, stencil_container_type, state_t, MeshType, ghost_container_type>;
    stencil_filler_t FillStencilX(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				   state, m_meshObj.get(), m_ghostLeft, m_ghostRight,
				   stencilVals, xAxis);

    stencil_filler_t FillStencilY(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				   state, m_meshObj.get(), m_ghostBack, m_ghostFront,
				   stencilVals, yAxis);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impl::ComputeDirectionalFluxBalanceFirstOrderJacobianOnBoundaryCell<
	  pda::impladvdiff2d::ComputeDirectionalFluxValuesAndJacobians<
	    pda::impl::ReconstructorFromStencil<
	      edge_rec_type, stencil_container_type>,
	    scalar_type, flux_type, flux_jac_type>,
	  dimensionality, MeshType, jacobian_type>,
      V_t, scalar_type
      >;

    functor_type funcx(V, m_meshObj.get().dxInv(),
		       /* end args for velo */
		       J, xAxis, m_meshObj.get(),
		       /* end args for jac */
		       m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR,
		       fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
		       /* end args for flux */
		       toReconstructionScheme(m_inviscidFluxRecEn), stencilVals,
		       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		       /* end args for reconstructor */
		       );

    functor_type funcy(V, m_meshObj.get().dyInv(),
		       /* end args for velo */
		       J, yAxis, m_meshObj.get(),
		       /* end args for jac */
		       m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF,
		       fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
		       /* end args for flux */
		       toReconstructionScheme(m_inviscidFluxRecEn), stencilVals,
		       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		       /* end args for reconstructor */
		       );

    const auto & graph      = m_meshObj.get().graph();
    constexpr auto two      = static_cast<scalar_type>(2);
    const auto & graphRows  = m_meshObj.get().graphRowsOfCellsNearBd();
    const auto dxInvSq	    = m_meshObj.get().dxInv()*m_meshObj.get().dxInv();
    const auto dyInvSq	    = m_meshObj.get().dyInv()*m_meshObj.get().dyInv();
    const auto diffDxInvSq  = m_burgers2d_diffusion*dxInvSq;
    const auto diffDyInvSq  = m_burgers2d_diffusion*dyInvSq;

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];
      const auto vIndex      = smPt*numDofPerCell;
      const auto uIndex	     = graph(smPt, 0)*numDofPerCell;
      const auto uIndexLeft  = graph(smPt, 1)*numDofPerCell;
      const auto uIndexFront = graph(smPt, 2)*numDofPerCell;
      const auto uIndexRight = graph(smPt, 3)*numDofPerCell;
      const auto uIndexBack  = graph(smPt, 4)*numDofPerCell;

      // x-direction inviscid contributions
      FillStencilX(smPt, it, numDofPerCell);
      if constexpr(std::is_same_v<BCFunctorsHolderType, impl::NoOperation<void>>){
        fillJacFactorsForCellBd(smPt, xAxis);
      }
      else {
        throw std::runtime_error("Custom BC Jacobian factors not implemented");
      }
      funcx(smPt, numDofPerCell, m_bcCellJacFactors);
      // u diffusion velocity contributions
      V(vIndex) += diffDxInvSq*( stencilVals(4) - two*stencilVals(2) + stencilVals(0) );
      V(vIndex) += diffDyInvSq*( stencilVals(4) - two*stencilVals(2) + stencilVals(0) );

      // y-direction inviscid contributions
      FillStencilY(smPt, it, numDofPerCell);
      if constexpr(std::is_same_v<BCFunctorsHolderType, impl::NoOperation<void>>){
	fillJacFactorsForCellBd(smPt, yAxis);
      }
      else {
        throw std::runtime_error("Custom BC Jacobian factors not implemented");
      }
      funcy(smPt, numDofPerCell, m_bcCellJacFactors);
      // y-direction diffusion velocity contributions
      V(vIndex+1) += diffDxInvSq*( stencilVals(5) - two*stencilVals(3) + stencilVals(1) );
      V(vIndex+1) += diffDyInvSq*( stencilVals(5) - two*stencilVals(3) + stencilVals(1) );

      // diffusion Jacobian contributions
      // self
      J.coeffRef(vIndex,   uIndex)   += -two*diffDxInvSq - two*diffDyInvSq;
      J.coeffRef(vIndex+1, uIndex+1) += -two*diffDxInvSq - two*diffDyInvSq;

      if (uIndexLeft > -1){
	J.coeffRef(vIndex,   uIndexLeft)   += diffDxInvSq;
        J.coeffRef(vIndex+1, uIndexLeft+1) += diffDxInvSq;
      }
      else{
        J.coeffRef(vIndex,   uIndex)   += -diffDxInvSq;
        J.coeffRef(vIndex+1, uIndex+1) += -diffDxInvSq;
      }

      if (uIndexFront > -1){
	J.coeffRef(vIndex,   uIndexFront)   += diffDyInvSq;
        J.coeffRef(vIndex+1, uIndexFront+1) += diffDyInvSq;
      }
      else{
        J.coeffRef(vIndex,   uIndex)   += -diffDyInvSq;
        J.coeffRef(vIndex+1, uIndex+1) += -diffDyInvSq;
      }

      if (uIndexRight > -1){
	J.coeffRef(vIndex,   uIndexRight)   += diffDxInvSq;
        J.coeffRef(vIndex+1, uIndexRight+1) += diffDxInvSq;
      }
      else{
	J.coeffRef(vIndex,   uIndex)   += -diffDxInvSq;
        J.coeffRef(vIndex+1, uIndex+1) += -diffDxInvSq;
      }

      if (uIndexBack > -1){
	J.coeffRef(vIndex,   uIndexBack)   += diffDyInvSq;
        J.coeffRef(vIndex+1, uIndexBack+1) += diffDyInvSq;
      }
      else{
	J.coeffRef(vIndex,   uIndex)   += -diffDyInvSq;
        J.coeffRef(vIndex+1, uIndex+1) += -diffDyInvSq;
      }
    }
  }

  template<class state_t, class V_t>
  void velocityOnlyNearBdCellsImpl(const state_t & state,
				   const scalar_type currentTime,
				   V_t & V,
				   flux_type & fluxL,
				   flux_type & fluxF,
				   flux_type & fluxR,
				   flux_type & fluxB,
				   edge_rec_type & uMinusHalfNeg,
				   edge_rec_type & uMinusHalfPos,
				   edge_rec_type & uPlusHalfNeg,
				   edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;

    const auto stencilSize = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    stencil_container_type stencilVals(numDofPerCell*stencilSize);

    // stencil filler needed because we are doing cells near boundaries
    using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, stencil_container_type, state_t, MeshType, ghost_container_type>;

    sfiller_t StencilFillerX(stencilSize, state, m_meshObj.get(), m_ghostLeft, m_ghostRight, stencilVals, xAxis);
    sfiller_t StencilFillerY(stencilSize, state, m_meshObj.get(), m_ghostBack, m_ghostFront, stencilVals, yAxis);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impladvdiff2d::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  scalar_type, flux_type>,
      V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.get().dxInv(),
		    /* end args for velo */
		    m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR,
		    /* end args for flux */
		    toReconstructionScheme(m_inviscidFluxRecEn), stencilVals,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.get().dyInv(),
		    /* end args for velo */
		    m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF,
		    /* end args for flux */
		    toReconstructionScheme(m_inviscidFluxRecEn), stencilVals,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		    /* end args for reconstructor */
		    );

    constexpr auto two      = static_cast<scalar_type>(2);
    const auto dxInvSq	    = m_meshObj.get().dxInv()*m_meshObj.get().dxInv();
    const auto dyInvSq	    = m_meshObj.get().dyInv()*m_meshObj.get().dyInv();
    const auto diffDxInvSq  = m_burgers2d_diffusion*dxInvSq;
    const auto diffDyInvSq  = m_burgers2d_diffusion*dyInvSq;

    const auto & rows   = m_meshObj.get().graphRowsOfCellsNearBd();
#if defined PRESSIODEMOAPPS_ENABLE_OPENMP && !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
#pragma omp for schedule(static)
#endif
    for (std::size_t it=0; it<rows.size(); ++it)
    {
      const auto smPt        = rows[it];
      const auto vIndex      = smPt*numDofPerCell;

      StencilFillerX(smPt, it, numDofPerCell);
      Fx(smPt, numDofPerCell);
      // *** add X contribution of diffusion ***
      if (stencilSize == 3){
	V(vIndex) += diffDxInvSq*( stencilVals(4) - two*stencilVals(2) + stencilVals(0) );
        V(vIndex) += diffDyInvSq*( stencilVals(4) - two*stencilVals(2) + stencilVals(0) );
      }
      else if (stencilSize == 5){
        throw std::runtime_error("Higher order not fixed yet");
	V(vIndex) += diffDxInvSq*( stencilVals(3) -two*stencilVals(2) +stencilVals(1) );
      }
      else if (stencilSize == 7){
        throw std::runtime_error("Higher order not fixed yet");
	V(vIndex) += diffDxInvSq*( stencilVals(4) -two*stencilVals(3) +stencilVals(2) );
      }

      StencilFillerY(smPt, it, numDofPerCell);
      Fy(smPt, numDofPerCell);
      // *** add Y contribution of diffusion ***
      if (stencilSize == 3){
	V(vIndex+1) += diffDxInvSq*( stencilVals(5) - two*stencilVals(3) + stencilVals(1) );
        V(vIndex+1) += diffDyInvSq*( stencilVals(5) - two*stencilVals(3) + stencilVals(1) );
      }
      else if (stencilSize == 5){
        throw std::runtime_error("Higher order not fixed yet");
	V(smPt) += diffDyInvSq*( stencilVals(3) -two*stencilVals(2) +stencilVals(1) );
      }
      else if (stencilSize == 7){
        throw std::runtime_error("Higher order not fixed yet");
	V(smPt) += diffDyInvSq*( stencilVals(4) -two*stencilVals(3) +stencilVals(2) );
      }
    }
  }


  template<class U_t, class V_t>
  void velocityOnlyInnerCellsImpl(const U_t & U,
				  const scalar_type /*currentTime*/,
				  V_t & V,
				  flux_type & fluxL,
				  flux_type & fluxF,
				  flux_type & fluxR,
				  flux_type & fluxB,
				  edge_rec_type & uMinusHalfNeg,
				  edge_rec_type & uMinusHalfPos,
				  edge_rec_type & uPlusHalfNeg,
				  edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impladvdiff2d::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorForDiscreteFunction<
	    dimensionality, MeshType, U_t, edge_rec_type>,
	  scalar_type, flux_type>,
      V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.get().dxInv(),
		    /* end args for velo */
		    m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR,
		    /* end args for flux */
		    xAxis, toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj.get(),
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.get().dyInv(),
		    /* end args for velo */
		    m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF,
		    /* end args for flux */
		    yAxis, toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj.get(),
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		    /* end args for reconstructor */
		    );

    const auto & graphRows = m_meshObj.get().graphRowsOfCellsAwayFromBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      Fx(smPt, numDofPerCell);
      Fy(smPt, numDofPerCell);

      // diffusion contribution
      addBurgersDiffusionToVelocityAndOptionalJacobianInnerCells(U, V, nullptr, smPt);
    }
  }

  void fillJacFactorsForCellBd(index_t graphRow, int axis) const
  {
    assert(axis <= 2);

    if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersOutflow) {
        // homogeneous Neumann
        m_bcCellJacFactors.fill(static_cast<scalar_type>(1));
    }
    else if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersPeriodic) {
        throw std::runtime_error("Should not be getting Jacobian factors for fully periodic problem.");
    }
    else {
        throw std::runtime_error("Invalid problem enumeration.");
    }
  }

  template<class U_t, class V_t>
  void addBurgersDiffusionToVelocityAndOptionalJacobianInnerCells(const U_t & U,
								   V_t & V,
								   jacobian_type * J,
								   index_t smPt) const
  {
    constexpr auto two  = static_cast<scalar_type>(2);
    const auto dxInvSq  = m_meshObj.get().dxInv()*m_meshObj.get().dxInv();
    const auto dyInvSq  = m_meshObj.get().dyInv()*m_meshObj.get().dyInv();
    const auto diffDxInvSq  = m_burgers2d_diffusion*dxInvSq;
    const auto diffDyInvSq  = m_burgers2d_diffusion*dyInvSq;

    const auto & graph     = m_meshObj.get().graph();
    const auto vIndex      = smPt*numDofPerCell;
    const auto uIndex      = graph(smPt, 0)*numDofPerCell;
    const auto uIndexLeft  = graph(smPt, 1)*numDofPerCell;
    const auto uIndexFront = graph(smPt, 2)*numDofPerCell;
    const auto uIndexRight = graph(smPt, 3)*numDofPerCell;
    const auto uIndexBack  = graph(smPt, 4)*numDofPerCell;

    V(vIndex) += diffDxInvSq*( U(uIndexRight) - two*U(uIndex) + U(uIndexLeft) );
    V(vIndex) += diffDyInvSq*( U(uIndexFront) - two*U(uIndex) + U(uIndexBack) );

    V(vIndex+1) += diffDxInvSq*( U(uIndexRight+1) - two*U(uIndex+1) + U(uIndexLeft+1) );
    V(vIndex+1) += diffDyInvSq*( U(uIndexFront+1) - two*U(uIndex+1) + U(uIndexBack+1) );

    if (J){
      J->coeffRef(vIndex, uIndex) += -two*diffDxInvSq - two*diffDyInvSq;
      J->coeffRef(vIndex, uIndexLeft)  += diffDxInvSq;
      J->coeffRef(vIndex, uIndexFront) += diffDyInvSq;
      J->coeffRef(vIndex, uIndexRight) += diffDxInvSq;
      J->coeffRef(vIndex, uIndexBack)  += diffDyInvSq;
      J->coeffRef(vIndex+1, uIndex+1) += -two*diffDxInvSq - two*diffDyInvSq;
      J->coeffRef(vIndex+1, uIndexLeft+1)  += diffDxInvSq;
      J->coeffRef(vIndex+1, uIndexFront+1) += diffDyInvSq;
      J->coeffRef(vIndex+1, uIndexRight+1) += diffDxInvSq;
      J->coeffRef(vIndex+1, uIndexBack+1)  += diffDyInvSq;
    }
  }

  void allocateGhosts()
  {
    const auto stencilSize    = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    const auto numGhostValues = numDofPerCell*((stencilSize-1)/2);

    const index_t s1 = m_meshObj.get().numCellsNearBd();
    ::pressiodemoapps::resize(m_ghostLeft, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostFront,s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight,s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostBack, s1, numGhostValues);
  }

protected:
  // common to all problems
  ::pressiodemoapps::AdvectionDiffusion2d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_inviscidFluxRecEn;
  ::pressiodemoapps::InviscidFluxScheme m_inviscidFluxSchemeEn;
  ::pressiodemoapps::ViscousFluxReconstruction m_viscousFluxRecEn;
//   BCFunctorsHolderType m_bcFuncsHolder = {};

  std::reference_wrapper<const MeshType> m_meshObj;
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  mutable ghost_container_type m_ghostLeft;
  mutable ghost_container_type m_ghostFront;
  mutable ghost_container_type m_ghostRight;
  mutable ghost_container_type m_ghostBack;

  std::array<scalar_type, 2> normalX_{1, 0};
  std::array<scalar_type, 2> normalY_{0, 1};

  // parameters specific to problems
  // will need to handle this better later
  scalar_type m_burgers2d_icPulse = {};
  scalar_type m_burgers2d_icSpread = {};
  scalar_type m_burgers2d_diffusion = {};
  scalar_type m_burgers2d_x0 = {};
  scalar_type m_burgers2d_y0 = {};

  mutable std::array<scalar_type, numDofPerCell> m_bcCellJacFactors;
};

template<class T1, class T2> constexpr int EigenApp<T1,T2>::numDofPerCell;
template<class T1, class T2> constexpr int EigenApp<T1,T2>::dimensionality;

}}//end namespace
#endif
