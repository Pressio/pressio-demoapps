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

#ifndef PRESSIODEMOAPPS_ADVECTION_DIFFUSION_REACTION_2D_PROB_CLASS_HPP_
#define PRESSIODEMOAPPS_ADVECTION_DIFFUSION_REACTION_2D_PROB_CLASS_HPP_

#include "advection_diffusion_reaction_2d_initial_condition.hpp"
#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "advection_diffusion_reaction_2d_flux_mixin.hpp"
#include "advection_diffusion_reaction_2d_ghost_filler_problemA.hpp"
#include "mixin_directional_flux_balance.hpp"
#include "mixin_directional_flux_balance_jacobian.hpp"
#include "Eigen/Sparse"

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#include <omp.h>
#endif

namespace pressiodemoapps{
namespace impladvdiffreac2d{

// tags are used inside he public create function: create_problem_...()
// in the file ../advection_diffusion.hpp
// to dispatch to the proper problem
// so add new ones if a new problem is added
struct TagProblemA{};

template<class MeshType>
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

private:
  static constexpr int dimensionality{2};

  using ghost_container_type      = Eigen::Matrix<scalar_type,
						  Eigen::Dynamic,
						  Eigen::Dynamic,
						  Eigen::RowMajor>;
  using stencil_container_type    = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using flux_type	          = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using edge_rec_type	          = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using flux_jac_type             = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;
  using reconstruction_gradient_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;

public:

  //
  // constructor for ProblemA
  //
  template<class SourceT>
  EigenApp(TagProblemA /*tag*/,
	   const MeshType & meshObj,
	   ::pressiodemoapps::InviscidFluxReconstruction inviscidFluxRecEn,
	   ::pressiodemoapps::InviscidFluxScheme invFluxSchemeEn,
	   ::pressiodemoapps::ViscousFluxReconstruction visFluxRecEn,
	   scalar_type ux,
	   scalar_type uy,
	   scalar_type diffusion,
	   scalar_type sigma_reaction,
	   SourceT sf)
    : m_numDofPerCell(1),
      m_probEn(::pressiodemoapps::AdvectionDiffusionReaction2d::ProblemA),
      m_inviscidFluxRecEn(inviscidFluxRecEn),
      m_inviscidFluxSchemeEn(invFluxSchemeEn),
      m_viscousFluxRecEn(visFluxRecEn),
      m_meshObj(meshObj),
      m_probA_ux(ux),
      m_probA_uy(uy),
      m_probA_diff(diffusion),
      m_probA_sigma_reaction(sigma_reaction)
  {
    m_probA_sourceFunctor  = sf;
    m_numDofStencilMesh = m_meshObj.stencilMeshSize() * m_numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * m_numDofPerCell;
    allocateGhosts();
  }

  state_type initialCondition() const{
    state_type initialState(m_numDofStencilMesh);
    if (m_probEn == ::pressiodemoapps::AdvectionDiffusionReaction2d::ProblemA){
      problemAIC(initialState, m_meshObj);
    }

    return initialState;
  }

protected:
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
      edge_rec_type uMinusHalfNeg(m_numDofPerCell);
      edge_rec_type uMinusHalfPos(m_numDofPerCell);
      edge_rec_type uPlusHalfNeg(m_numDofPerCell);
      edge_rec_type uPlusHalfPos(m_numDofPerCell);
      // fluxes
      flux_type fluxL(m_numDofPerCell);
      flux_type fluxF(m_numDofPerCell);;
      flux_type fluxR(m_numDofPerCell);
      flux_type fluxB(m_numDofPerCell);

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

      fillGhosts(U);

      if (J){
	velocityAndJacobianImpl(U, currentTime, V, *J,
				fluxL, fluxF, fluxR, fluxB,
				uMinusHalfNeg, uMinusHalfPos,
				uPlusHalfNeg,  uPlusHalfPos);
      }
      else{
	velocityOnlyNearBdCellsImpl(U, currentTime, V,
				    fluxL, fluxF, fluxR, fluxB,
				    uMinusHalfNeg, uMinusHalfPos,
				    uPlusHalfNeg,  uPlusHalfPos);

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
  void fillGhosts(const U_t & U) const
  {
    if (m_probEn != ::pressiodemoapps::AdvectionDiffusionReaction2d::ProblemA){
      throw std::runtime_error("fillGhosts not imple for this problem");
    }

    const auto stencilSize = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    using ghost_filler_t  = ::pressiodemoapps::impladvdiffreac2d::GhostFillerProblemA<
      U_t, MeshType, ghost_container_type>;
    ghost_filler_t ghF(stencilSize, U, m_meshObj,
		       m_ghostLeft, m_ghostFront,
		       m_ghostRight, m_ghostBack);
    const auto & rowsBd = m_meshObj.graphRowsOfCellsNearBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(rowsBd.size()) it=0; it<rowsBd.size(); ++it){
      ghF(rowsBd[it], it);
    }
  }

  template<class Tr>
  void initializeJacobianForInnerCells(std::vector<Tr> & trList)
  {
    // for inner cells, the Jacobian is exact and depends on
    // the scheme wanted by the user, no special treatment needed

    const scalar_type zero = 0;
    const auto & graph = m_meshObj.graph();
    // only grab the graph rows for INNER cells (i.e. AWAY from boundaries)
    const auto & targetGraphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (std::size_t it=0; it<targetGraphRows.size(); ++it)
      {
	const auto smPt = targetGraphRows[it];

	// find out which row in the jacobian we are dealing with
	const auto jacRowOfCurrCellFirstDof = smPt*m_numDofPerCell;

	// initialize jacobian block entries wrt current cell's dofs
	const auto jacColOfCurrCellRho = graph(smPt, 0)*m_numDofPerCell;
	for (int k=0; k<m_numDofPerCell; ++k){
	  for (int j=0; j<m_numDofPerCell; ++j){
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
	  const auto colInd = graph(smPt, i)*m_numDofPerCell;
	  for (int k=0; k<m_numDofPerCell; ++k){
	    for (int j=0; j<m_numDofPerCell; ++j){
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
    const auto & graph = m_meshObj.graph();
    const auto & targetGraphRows = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<targetGraphRows.size(); ++it)
      {
	const auto smPt = targetGraphRows[it];
	const auto jacRowOfCurrCellRho = smPt*m_numDofPerCell;
	const auto jacColOfCurrCellRho = graph(smPt, 0)*m_numDofPerCell;

	// wrt current cell's dofs
	for (int k=0; k<m_numDofPerCell; ++k){
	  for (int j=0; j<m_numDofPerCell; ++j){
	    trList.push_back( Tr(jacRowOfCurrCellRho+k, jacColOfCurrCellRho+j, zero) );
	  }
	}

	// wrt neighbors
	// for near-bd, we only do first-order Jacobian for now
	for (int i=1; i<=4; ++i){
	  const auto nID = graph(smPt, i);
	  if (nID != -1){
	    const auto ci = nID*m_numDofPerCell;
	    for (int k=0; k<m_numDofPerCell; ++k){
	      for (int j=0; j<m_numDofPerCell; ++j){
		trList.push_back( Tr(jacRowOfCurrCellRho+k, ci+j, zero) );
	      }
	    }
	  }
	}
      }
  }


  template<class state_t, class V_t>
  void velocityAndJacobianImpl(const state_t & state,
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
    flux_jac_type fluxJacLNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacLPos(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacFNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacFPos(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacRNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacRPos(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacBNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacBPos(m_numDofPerCell, m_numDofPerCell);

    int nonZerosCountBeforeComputing = J.nonZeros();
    (void) nonZerosCountBeforeComputing;

    // near boundary I have be careful because
    // the jacobian can only be first order for now
    if (m_inviscidFluxRecEn == InviscidFluxReconstruction::FirstOrder){
      velocityAndJacNearBDCellsImplFirstOrder(state, currentTime, V, J,
					      fluxL, fluxF, fluxR, fluxB,
					      fluxJacLNeg, fluxJacLPos,
					      fluxJacFNeg, fluxJacFPos,
					      fluxJacRNeg, fluxJacRPos,
					      fluxJacBNeg, fluxJacBPos,
					      uMinusHalfNeg, uMinusHalfPos,
					      uPlusHalfNeg,  uPlusHalfPos);
    }
    else{
      velocityAndJacNearBDCellsImplDifferentSchemeNotFirstOrderInviscid
	(state, currentTime, V, J,
	 fluxL, fluxF, fluxR, fluxB,
	 fluxJacLNeg, fluxJacLPos,
	 fluxJacFNeg, fluxJacFPos,
	 fluxJacRNeg, fluxJacRPos,
	 fluxJacBNeg, fluxJacBPos,
	 uMinusHalfNeg, uMinusHalfPos,
	 uPlusHalfNeg,  uPlusHalfPos);
    }

    velocityAndJacInnerCellsImpl(state, currentTime, V, J,
				 fluxL, fluxF, fluxR, fluxB,
				 fluxJacLNeg, fluxJacLPos,
				 fluxJacFNeg, fluxJacFPos,
				 fluxJacRNeg, fluxJacRPos,
				 fluxJacBNeg, fluxJacBPos,
				 uMinusHalfNeg, uMinusHalfPos,
				 uPlusHalfNeg,  uPlusHalfPos);

    assert(J.nonZeros() == nonZerosCountBeforeComputing);
  }

  template<class state_t, class V_t>
  void velocityAndJacInnerCellsImpl(const state_t & state,
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
    // for inner cells, velocity and Jacobian
    // are both computed according to the order selected by the user
    // because for inner cells we support also Jacobians for Weno

    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;

    // allocate gradients of reconstructed states
    // the size depends on the scheme selected
    const auto stencilSize = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    reconstruction_gradient_t gradLNeg(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradLPos(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradFNeg(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradFPos(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRNeg(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRPos(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradBNeg(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradBPos(m_numDofPerCell, stencilSize-1);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impl::ComputeDirectionalFluxBalanceJacobianOnInteriorCell<
	  pda::impladvdiffreac2d::ComputeDirectionalFluxValuesAndJacobians<
	    pda::impl::ReconstructorForDiscreteFunction<
	      dimensionality, MeshType, state_t, edge_rec_type, reconstruction_gradient_t>,
	    scalar_type, flux_type, flux_jac_type>,
	  dimensionality, MeshType, jacobian_type>,
      V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.dxInv(),
		    /* end args for velo */
		    J, xAxis, m_meshObj,
		    /* end args for jac */
		    m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR,
		    fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos, m_probA_ux,
		    /* end args for flux */
		    xAxis, toReconstructionScheme(m_inviscidFluxRecEn), state, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		    gradLNeg, gradLPos, gradRNeg, gradRPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.dyInv(),
		    /* end args for velo */
		    J, yAxis, m_meshObj,
		    /* end args for jac */
		    m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF,
		    fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos, m_probA_uy,
		    /* end args for flux */
		    yAxis, toReconstructionScheme(m_inviscidFluxRecEn), state, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		    gradBNeg, gradBPos, gradFNeg, gradFPos
		    /* end args for reconstructor */
		    );

    const auto & x         = m_meshObj.viewX();
    const auto & y         = m_meshObj.viewY();
    const auto & graph     = m_meshObj.graph();
    const auto & graphRows = m_meshObj.graphRowsOfCellsAwayFromBd();

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    scalar_type sourceVal = {};
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];
      const auto uIndex	= graph(smPt, 0);

      Fx(smPt, m_numDofPerCell);
      Fy(smPt, m_numDofPerCell);

      addDiffusionToVelocityAndOptionalJacobianInnerCells(state, V, &J, smPt);

      // source term: this will overwrite sourceVal
      m_probA_sourceFunctor(x(uIndex), y(uIndex), currentTime, sourceVal);
      V(smPt) += sourceVal;

      // reaction contribution (note the minus sign)
      V(smPt) -= m_probA_sigma_reaction * state(uIndex);
      const auto vIndex = smPt*m_numDofPerCell;
      J.coeffRef(vIndex, uIndex) -= m_probA_sigma_reaction;
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

    stencil_container_type stencilValsForV(m_numDofPerCell*stencilSizeForV);

    stencil_filler_t FillStencilVeloX(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				      state, m_meshObj, m_ghostLeft, m_ghostRight,
				      stencilValsForV, xAxis);
    stencil_filler_t FillStencilVeloY(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				      state, m_meshObj, m_ghostBack, m_ghostFront,
				      stencilValsForV, yAxis);

    using velo_functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impladvdiffreac2d::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  scalar_type, flux_type>,
      V_t, scalar_type
      >;

    velo_functor_type funcVeloX(V, m_meshObj.dxInv(),
				/* end args for velo */
				m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR, m_probA_ux,
				/* end args for flux */
				toReconstructionScheme(m_inviscidFluxRecEn), stencilValsForV,
				uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
				/* end args for reconstructor */
				);

    velo_functor_type funcVeloY(V, m_meshObj.dyInv(),
				/* end args for velo */
				m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF, m_probA_uy,
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
    stencil_container_type stencilValsForJ(m_numDofPerCell*stencilSizeForJ);
    stencil_filler_t FillStencilJacX(stencilSizeForJ,
				     state, m_meshObj, m_ghostLeft, m_ghostRight,
				     stencilValsForJ, xAxis);
    stencil_filler_t FillStencilJacY(stencilSizeForJ,
				     state, m_meshObj, m_ghostBack, m_ghostFront,
				     stencilValsForJ, yAxis);

    using jac_functor_type =
      pda::impl::ComputeDirectionalFluxBalanceFirstOrderJacobianOnBoundaryCell<
	pda::impladvdiffreac2d::ComputeDirectionalFluxJacobians<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  scalar_type, flux_jac_type>,
      dimensionality, MeshType, jacobian_type
      >;

    jac_functor_type funcJacX(J, xAxis, m_meshObj,
			      /* end args for jac */
			      m_probEn, m_inviscidFluxSchemeEn, normalX_,
			      fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos, m_probA_ux,
			      /* end args for flux */
			      toReconstructionScheme(firstOrderRec), stencilValsForJ,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			      /* end args for reconstructor */
			      );

    jac_functor_type funcJacY(J, yAxis, m_meshObj,
			      /* end args for jac */
			      m_probEn, m_inviscidFluxSchemeEn, normalY_,
			      fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos, m_probA_uy,
			      /* end args for flux */
			      toReconstructionScheme(firstOrderRec), stencilValsForJ,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			      /* end args for reconstructor */
			      );

    std::array<scalar_type, 1> bdCellJacFactorsX{-1};
    std::array<scalar_type, 1> bdCellJacFactorsY{-1};

    // ************
    // loop
    // ************
    const auto & x	   = m_meshObj.viewX();
    const auto & y         = m_meshObj.viewY();
    const auto & graph     = m_meshObj.graph();
    const auto & graphRows = m_meshObj.graphRowsOfCellsNearBd();
    const auto dxInvSq	   = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto dyInvSq	   = m_meshObj.dyInv()*m_meshObj.dyInv();
    const auto diffDxInvSq = m_probA_diff*dxInvSq;
    const auto diffDyInvSq = m_probA_diff*dyInvSq;
    constexpr auto two      = static_cast<scalar_type>(2);

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    scalar_type sourceVal = {};
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];
      const auto uIndex	= graph(smPt, 0);
      const auto uIndexLeft  = graph(smPt, 1);
      const auto uIndexFront = graph(smPt, 2);
      const auto uIndexRight = graph(smPt, 3);
      const auto uIndexBack  = graph(smPt, 4);

      FillStencilVeloX(smPt, it, m_numDofPerCell);
      funcVeloX(smPt, m_numDofPerCell);
      FillStencilJacX(smPt, it, m_numDofPerCell);
      funcJacX(smPt, m_numDofPerCell, bdCellJacFactorsX);
      // case == 3 cannot happen here, see comment at near top of this func
      if (stencilSizeForV == 5){
	V(smPt) += diffDxInvSq*( stencilValsForV(3) -two*stencilValsForV(2) +stencilValsForV(1) );
      }
      else if (stencilSizeForV == 7){
	V(smPt) += diffDxInvSq*( stencilValsForV(4) -two*stencilValsForV(3) +stencilValsForV(2) );
      }

      FillStencilVeloY(smPt, it, m_numDofPerCell);
      funcVeloY(smPt, m_numDofPerCell);
      FillStencilJacY(smPt, it, m_numDofPerCell);
      funcJacY(smPt, m_numDofPerCell, bdCellJacFactorsY);
      if (stencilSizeForV == 5){
	V(smPt) += diffDyInvSq*( stencilValsForV(3) -two*stencilValsForV(2) +stencilValsForV(1) );
      }
      else if (stencilSizeForV == 7){
	V(smPt) += diffDyInvSq*( stencilValsForV(4) -two*stencilValsForV(3) +stencilValsForV(2) );
      }

      // source term: this will overwrite sourceVal
      m_probA_sourceFunctor(x(uIndex), y(uIndex), currentTime,  sourceVal);
      V(smPt) += sourceVal;

      // add reaction contribution
      V(smPt) -= m_probA_sigma_reaction * state(uIndex);

      auto selfValue = -two*diffDxInvSq -two*diffDyInvSq - m_probA_sigma_reaction;
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
    stencil_container_type stencilVals(m_numDofPerCell*stencilSize);

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, stencil_container_type, state_t, MeshType, ghost_container_type>;
    stencil_filler_t FillStencilX(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				   state, m_meshObj, m_ghostLeft, m_ghostRight,
				   stencilVals, xAxis);

    stencil_filler_t FillStencilY(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				   state, m_meshObj, m_ghostBack, m_ghostFront,
				   stencilVals, yAxis);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impl::ComputeDirectionalFluxBalanceFirstOrderJacobianOnBoundaryCell<
	  pda::impladvdiffreac2d::ComputeDirectionalFluxValuesAndJacobians<
	    pda::impl::ReconstructorFromStencil<
	      edge_rec_type, stencil_container_type>,
	    scalar_type, flux_type, flux_jac_type>,
	  dimensionality, MeshType, jacobian_type>,
      V_t, scalar_type
      >;

    functor_type funcx(V, m_meshObj.dxInv(),
		       /* end args for velo */
		       J, xAxis, m_meshObj,
		       /* end args for jac */
		       m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR,
		       fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos, m_probA_ux,
		       /* end args for flux */
		       toReconstructionScheme(m_inviscidFluxRecEn), stencilVals,
		       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		       /* end args for reconstructor */
		       );

    functor_type funcy(V, m_meshObj.dyInv(),
		       /* end args for velo */
		       J, yAxis, m_meshObj,
		       /* end args for jac */
		       m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF,
		       fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos, m_probA_uy,
		       /* end args for flux */
		       toReconstructionScheme(m_inviscidFluxRecEn), stencilVals,
		       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		       /* end args for reconstructor */
		       );

    std::array<scalar_type, 1> bdCellJacFactorsX;
    std::array<scalar_type, 1> bdCellJacFactorsY;
    bdCellJacFactorsX.fill(static_cast<scalar_type>(-1));
    bdCellJacFactorsY.fill(static_cast<scalar_type>(-1));

    const auto & x	    = m_meshObj.viewX();
    const auto & y          = m_meshObj.viewY();
    const auto & graph      = m_meshObj.graph();
    constexpr auto two      = static_cast<scalar_type>(2);
    const auto & graphRows  = m_meshObj.graphRowsOfCellsNearBd();
    const auto dxInvSq	    = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto dyInvSq	    = m_meshObj.dyInv()*m_meshObj.dyInv();
    const auto diffDxInvSq  = m_probA_diff*dxInvSq;
    const auto diffDyInvSq  = m_probA_diff*dyInvSq;

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    scalar_type sourceVal = {};
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];
      const auto uIndex	= graph(smPt, 0);
      const auto uIndexLeft  = graph(smPt, 1);
      const auto uIndexFront = graph(smPt, 2);
      const auto uIndexRight = graph(smPt, 3);
      const auto uIndexBack  = graph(smPt, 4);

      FillStencilX(smPt, it, m_numDofPerCell);
      funcx(smPt, m_numDofPerCell, bdCellJacFactorsX);
      V(smPt) += diffDxInvSq*( stencilVals(2) -two*stencilVals(1) +stencilVals(0) );

      FillStencilY(smPt, it, m_numDofPerCell);
      funcy(smPt, m_numDofPerCell, bdCellJacFactorsY);
      V(smPt) += diffDyInvSq*( stencilVals(2) -two*stencilVals(1) +stencilVals(0) );

      // source term: this will overwrite sourceVal
      m_probA_sourceFunctor(x(uIndex), y(uIndex), currentTime,  sourceVal);
      V(smPt) += sourceVal;

      // add reaction contribution
      V(smPt) -= m_probA_sigma_reaction * state(uIndex);

      auto selfValue = -two*diffDxInvSq -two*diffDyInvSq - m_probA_sigma_reaction;
      if (uIndexLeft != -1){
	J.coeffRef(smPt, uIndexLeft) += diffDxInvSq;
      }else{
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
    stencil_container_type stencilVals(stencilSize);

    // stencil filler needed because we are doing cells near boundaries
    using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, stencil_container_type, state_t, MeshType, ghost_container_type>;

    sfiller_t StencilFillerX(stencilSize, state, m_meshObj, m_ghostLeft, m_ghostRight, stencilVals, xAxis);
    sfiller_t StencilFillerY(stencilSize, state, m_meshObj, m_ghostBack, m_ghostFront, stencilVals, yAxis);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impladvdiffreac2d::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  scalar_type, flux_type>,
      V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.dxInv(),
		    /* end args for velo */
		    m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR, m_probA_ux,
		    /* end args for flux */
		    toReconstructionScheme(m_inviscidFluxRecEn), stencilVals,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.dyInv(),
		    /* end args for velo */
		    m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF, m_probA_uy,
		    /* end args for flux */
		    toReconstructionScheme(m_inviscidFluxRecEn), stencilVals,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		    /* end args for reconstructor */
		    );

    const auto & x	    = m_meshObj.viewX();
    const auto & y          = m_meshObj.viewY();
    const auto & graph      = m_meshObj.graph();
    constexpr auto two      = static_cast<scalar_type>(2);
    const auto dxInvSq	    = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto dyInvSq	    = m_meshObj.dyInv()*m_meshObj.dyInv();
    const auto diffDxInvSq  = m_probA_diff*dxInvSq;
    const auto diffDyInvSq  = m_probA_diff*dyInvSq;

    const auto & rows   = m_meshObj.graphRowsOfCellsNearBd();
#if defined PRESSIODEMOAPPS_ENABLE_OPENMP && !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
#pragma omp for schedule(static)
#endif
    scalar_type sourceVal = {};
    for (std::size_t it=0; it<rows.size(); ++it)
    {
      const auto smPt        = rows[it];
      const auto uIndex	     = graph(smPt, 0);

      // *** add X contribution of diffusion ***
      StencilFillerX(smPt, it, m_numDofPerCell);
      Fx(smPt, m_numDofPerCell);
      if (stencilSize == 3){
	V(smPt) += diffDxInvSq*( stencilVals(2) -two*stencilVals(1) +stencilVals(0) );
      }
      else if (stencilSize == 5){
	V(smPt) += diffDxInvSq*( stencilVals(3) -two*stencilVals(2) +stencilVals(1) );
      }
      else if (stencilSize == 7){
	V(smPt) += diffDxInvSq*( stencilVals(4) -two*stencilVals(3) +stencilVals(2) );
      }

      // *** add Y contribution of diffusion ***
      StencilFillerY(smPt, it, m_numDofPerCell);
      Fy(smPt, m_numDofPerCell);
      if (stencilSize == 3){
	V(smPt) += diffDyInvSq*( stencilVals(2) -two*stencilVals(1) +stencilVals(0) );
      }
      else if (stencilSize == 5){
	V(smPt) += diffDyInvSq*( stencilVals(3) -two*stencilVals(2) +stencilVals(1) );
      }
      else if (stencilSize == 7){
	V(smPt) += diffDyInvSq*( stencilVals(4) -two*stencilVals(3) +stencilVals(2) );
      }

      // source term: this will overwrite sourceVal
      m_probA_sourceFunctor(x(uIndex), y(uIndex), currentTime,  sourceVal);
      V(smPt) += sourceVal;

      // add reaction contribution
      V(smPt) -= m_probA_sigma_reaction * state(uIndex);
    }
  }


  template<class state_t, class V_t>
  void velocityOnlyInnerCellsImpl(const state_t & state,
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

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impladvdiffreac2d::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorForDiscreteFunction<
	    dimensionality, MeshType, state_t, edge_rec_type>,
	  scalar_type, flux_type>,
      V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.dxInv(),
		    /* end args for velo */
		    m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR, m_probA_ux,
		    /* end args for flux */
		    xAxis, toReconstructionScheme(m_inviscidFluxRecEn), state, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.dyInv(),
		    /* end args for velo */
		    m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF, m_probA_uy,
		    /* end args for flux */
		    yAxis, toReconstructionScheme(m_inviscidFluxRecEn), state, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		    /* end args for reconstructor */
		    );

    const auto & x         = m_meshObj.viewX();
    const auto & y         = m_meshObj.viewY();
    const auto & graph     = m_meshObj.graph();
    const auto & graphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    scalar_type sourceVal = {};
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];
      const auto uIndex	= graph(smPt, 0);

      Fx(smPt, m_numDofPerCell);
      Fy(smPt, m_numDofPerCell);

      // diffusion contribution
      addDiffusionToVelocityAndOptionalJacobianInnerCells(state, V, nullptr, smPt);

      // source term: this will overwrite sourceVal
      m_probA_sourceFunctor(x(uIndex), y(uIndex), currentTime, sourceVal);
      V(smPt) += sourceVal;

      // reaction contribution
      V(smPt) -= m_probA_sigma_reaction * state(uIndex);
    }
  }

  template<class state_t, class V_t>
  void addDiffusionToVelocityAndOptionalJacobianInnerCells(const state_t & state,
							   V_t & V,
							   jacobian_type * J,
							   index_t smPt) const
  {
    constexpr auto two     = static_cast<scalar_type>(2);
    const auto dxInvSq     = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto dyInvSq     = m_meshObj.dyInv()*m_meshObj.dyInv();
    const auto diffDxInvSq = m_probA_diff*dxInvSq;
    const auto diffDyInvSq = m_probA_diff*dyInvSq;

    const auto & graph     = m_meshObj.graph();
    const auto vIndex      = smPt*m_numDofPerCell;
    const auto uIndex      = graph(smPt, 0)*m_numDofPerCell;
    const auto uIndexLeft  = graph(smPt, 1)*m_numDofPerCell;
    const auto uIndexFront = graph(smPt, 2)*m_numDofPerCell;
    const auto uIndexRight = graph(smPt, 3)*m_numDofPerCell;
    const auto uIndexBack  = graph(smPt, 4)*m_numDofPerCell;
    V(vIndex)   += diffDxInvSq*( state(uIndexRight) - two*state(uIndex) + state(uIndexLeft) );
    V(vIndex)   += diffDyInvSq*( state(uIndexFront) - two*state(uIndex) + state(uIndexBack) );

    if (J){
      J->coeffRef(vIndex, uIndex)          += -two*diffDxInvSq - two*diffDyInvSq;
      J->coeffRef(vIndex, uIndexLeft)      += diffDxInvSq;
      J->coeffRef(vIndex, uIndexFront)     += diffDyInvSq;
      J->coeffRef(vIndex, uIndexRight)     += diffDxInvSq;
      J->coeffRef(vIndex, uIndexBack)      += diffDyInvSq;
    }
  }

  void allocateGhosts()
  {
    const auto stencilSize    = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    const auto numGhostValues = m_numDofPerCell*((stencilSize-1)/2);

    const index_t s1 = m_meshObj.numCellsBd();
    ::pressiodemoapps::resize(m_ghostLeft, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostFront,s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight,s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostBack, s1, numGhostValues);
  }

protected:
  // common to all problems
  int m_numDofPerCell = {};
  ::pressiodemoapps::AdvectionDiffusionReaction2d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_inviscidFluxRecEn;
  ::pressiodemoapps::InviscidFluxScheme m_inviscidFluxSchemeEn;
  ::pressiodemoapps::ViscousFluxReconstruction m_viscousFluxRecEn;

  const MeshType & m_meshObj;
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  mutable ghost_container_type m_ghostLeft;
  mutable ghost_container_type m_ghostFront;
  mutable ghost_container_type m_ghostRight;
  mutable ghost_container_type m_ghostBack;

  const std::array<scalar_type, 2> normalX_{1, 0};
  const std::array<scalar_type, 2> normalY_{0, 1};

  // parameters specific to problems
  // will need to handle this better later
  scalar_type m_probA_ux = {};
  scalar_type m_probA_uy = {};
  scalar_type m_probA_diff = {};
  scalar_type m_probA_sigma_reaction = {};
  std::function<void(const scalar_type & /*x*/,
		     const scalar_type & /*y*/,
		     const scalar_type & /*time*/,
		     scalar_type &)> m_probA_sourceFunctor;

};

template<class MeshType> constexpr int EigenApp<MeshType>::dimensionality;

}}//end namespace
#endif
