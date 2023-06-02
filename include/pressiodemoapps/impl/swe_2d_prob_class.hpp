/*
//@HEADER
// ************************************************************************
//
// swe_2d_prob_class.hpp
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

#ifndef PRESSIODEMOAPPS_SWE2D_IMPL_HPP_
#define PRESSIODEMOAPPS_SWE2D_IMPL_HPP_

#include "swe_rusanov_flux_values_function.hpp"
#include "swe_rusanov_flux_jacobian_function.hpp"
#include "swe_2d_initial_condition.hpp"
#include "swe_2d_ghost_filler_inviscid_wall.hpp"
#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "swe_2d_flux_mixin.hpp"
#include "mixin_directional_flux_balance.hpp"
#include "mixin_directional_flux_balance_jacobian.hpp"
#include "Eigen/Sparse"

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#include <omp.h>
#endif

namespace pressiodemoapps{
namespace implswe2d{

// tags are used inside he public create function: create_problem_...()
// in the file ../advection_diffusion.hpp
// to dispatch to the proper problem
// so add new ones if a new problem is added
struct TagProblemSlipWall{};


template<class MeshType>
class EigenApp
{

public:
  using index_t	      = typename MeshType::index_t;
  using scalar_type   = typename MeshType::scalar_t;
  using state_type    = Eigen::Matrix<scalar_type,Eigen::Dynamic,1>;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, index_t>;
  using graph_type    = typename MeshType::graph_t;

private:
  static constexpr int dimensionality{2};
  static constexpr int numDofPerCell{3};

  using ghost_container_type      = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using stencil_container_type    = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using flux_type	          = Eigen::Matrix<scalar_type, numDofPerCell,  1>;
  using edge_rec_type	          = Eigen::Matrix<scalar_type, numDofPerCell,  1>;
  using flux_jac_type             = Eigen::Matrix<scalar_type, numDofPerCell, numDofPerCell>;
  using reconstruction_gradient_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  EigenApp(TagProblemSlipWall /*tag*/,
	   const MeshType & meshObj,
	   ::pressiodemoapps::InviscidFluxReconstruction recEn,
	   ::pressiodemoapps::InviscidFluxScheme fluxEnum,
	   scalar_type gravity,
	   scalar_type coriolis,
	   scalar_type initialPulseMag)
    : m_probEn(::pressiodemoapps::Swe2d::SlipWall),
      m_recEn(recEn),
      m_fluxEn(fluxEnum),
      m_meshObj(meshObj),
      m_gravity(gravity),
      m_coriolis(coriolis),
      m_initialPulseMagnitude(initialPulseMag)
  {
    m_numDofStencilMesh = m_meshObj.stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * numDofPerCell;
    allocateGhosts();
  }

  scalar_type gravity()         const{ return m_gravity; }
  scalar_type coriolis()        const{ return m_coriolis; }

  state_type initialCondition() const
  {
    state_type initialState(m_numDofStencilMesh);
    switch(m_probEn){
    case ::pressiodemoapps::Swe2d::SlipWall:{
      GaussianPulse(initialState, m_meshObj, m_initialPulseMagnitude);
    }
    };

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

  // note that here we MUST use a template because when doing
  // bindings, this gets deduced to be a Ref
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

    // these are private variables for each thread
    // when openmp is enabled
    // reconstructions values
    edge_rec_type uMinusHalfNeg, uMinusHalfPos;
    edge_rec_type uPlusHalfNeg,  uPlusHalfPos;
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
  template<class Tr>
  void initializeJacobianForInnerCells(std::vector<Tr> & trList)
  {
    // for inner cells, the Jacobian is of the same scheme
    // wanted by the user, no special treatment is needed

    const scalar_type zero = 0;
    const auto & graph = m_meshObj.graph();
    const auto & targetGraphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (std::size_t it=0; it<targetGraphRows.size(); ++it)
      {
	const auto smPt = targetGraphRows[it];
	const auto jacRowOfCurrCellRho = smPt*numDofPerCell;

	// entries wrt current cell's dofs
	const auto jacColOfCurrCellRho = graph(smPt, 0)*numDofPerCell;
	for (int k=0; k<numDofPerCell; ++k){
	  for (int j=0; j<numDofPerCell; ++j){
	    trList.push_back( Tr(jacRowOfCurrCellRho+k, jacColOfCurrCellRho+j, zero) );
	  }
	}

	// wrt neighbors: this depends on the advection scheme
	const auto numNeighbors =
	  (m_recEn == InviscidFluxReconstruction::FirstOrder) ? 4
	    : (m_recEn == InviscidFluxReconstruction::Weno3) ? 8 : 12;

	for (int i=1; i<=numNeighbors; ++i){
	  const auto ci = graph(smPt, i)*numDofPerCell;
	  for (int k=0; k<numDofPerCell; ++k){
	    for (int j=0; j<numDofPerCell; ++j){
	      trList.push_back( Tr(jacRowOfCurrCellRho+k, ci+j, zero) );
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

  template<class U_t>
  void fillGhosts(const U_t & U) const
  {
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    if (m_probEn == ::pressiodemoapps::Swe2d::SlipWall)
    {
      using ghost_filler_t  = ::pressiodemoapps::implswe2d::InviscidWallFiller<
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

    else{
      // no op
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
    flux_jac_type fluxJacLNeg, fluxJacLPos;
    flux_jac_type fluxJacFNeg, fluxJacFPos;
    flux_jac_type fluxJacRNeg, fluxJacRPos;
    flux_jac_type fluxJacBNeg, fluxJacBPos;

    int nonZerosCountBeforeComputing = J.nonZeros();
    (void) nonZerosCountBeforeComputing;

    // near boundary I have be careful because
    // the jacobian can only be first order for now
    if (m_recEn == InviscidFluxReconstruction::FirstOrder){
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
      velocityAndJacNearBDCellsImplDifferentScheme(U, currentTime, V, J,
						   fluxL, fluxF, fluxR, fluxB,
						   fluxJacLNeg, fluxJacLPos,
						   fluxJacFNeg, fluxJacFPos,
						   fluxJacRNeg, fluxJacRPos,
						   fluxJacBNeg, fluxJacBPos,
						   uMinusHalfNeg, uMinusHalfPos,
						   uPlusHalfNeg,  uPlusHalfPos);
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
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
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
	  pda::implswe2d::ComputeDirectionalFluxValuesAndJacobians<
	    pda::impl::ReconstructorForDiscreteFunction<
	      dimensionality, MeshType, U_t, edge_rec_type, reconstruction_gradient_t>,
	    scalar_type, flux_type, flux_jac_type>,
	  dimensionality, MeshType, jacobian_type>,
      V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.dxInv(),
		    /* end args for velo */
		    J, xAxis, m_meshObj,
		    /* end args for jac */
		    m_fluxEn, normalX_, m_gravity, fluxL, fluxR,
		    fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
		    /* end args for flux */
		    xAxis, toReconstructionScheme(m_recEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		    gradLNeg, gradLPos, gradRNeg, gradRPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.dyInv(),
		    /* end args for velo */
		    J, yAxis, m_meshObj,
		    /* end args for jac */
		    m_fluxEn, normalY_, m_gravity, fluxB, fluxF,
		    fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
		    /* end args for flux */
		    yAxis, toReconstructionScheme(m_recEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		    gradBNeg, gradBPos, gradFNeg, gradFPos
		    /* end args for reconstructor */
		    );

    const auto & graphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];
      Fx(smPt, numDofPerCell);
      Fy(smPt, numDofPerCell);
      addForcingContributionToVelocityAndJacobian(U, V, J, smPt);
    }
  }

  template<class U_t, class V_t>
  void velocityAndJacNearBDCellsImplFirstOrder(const U_t & U,
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
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    assert(m_recEn == InviscidFluxReconstruction::FirstOrder);

    // if here, then the scheme for velocity matches
    // the one for Jacobian so we can use same functors

    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    stencil_container_type stencilVals(numDofPerCell*stencilSize);

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, stencil_container_type, U_t, MeshType, ghost_container_type>;
    stencil_filler_t FillStencilX(reconstructionTypeToStencilSize(m_recEn),
				   U, m_meshObj, m_ghostLeft, m_ghostRight,
				   stencilVals, xAxis);

    stencil_filler_t FillStencilY(reconstructionTypeToStencilSize(m_recEn),
				   U, m_meshObj, m_ghostBack, m_ghostFront,
				   stencilVals, yAxis);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impl::ComputeDirectionalFluxBalanceFirstOrderJacobianOnBoundaryCell<
	  pda::implswe2d::ComputeDirectionalFluxValuesAndJacobians<
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
		       m_fluxEn, normalX_, m_gravity, fluxL, fluxR,
		       fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
		       /* end args for flux */
		       toReconstructionScheme(m_recEn), stencilVals,
		       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		       /* end args for reconstructor */
		       );

    functor_type funcy(V, m_meshObj.dyInv(),
		       /* end args for velo */
		       J, yAxis, m_meshObj,
		       /* end args for jac */
		       m_fluxEn, normalY_, m_gravity, fluxB, fluxF,
		       fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
		       /* end args for flux */
		       toReconstructionScheme(m_recEn), stencilVals,
		       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		       /* end args for reconstructor */
		       );

    std::array<scalar_type, numDofPerCell> bdCellJacFactorsX;
    std::array<scalar_type, numDofPerCell> bdCellJacFactorsY;
    bdCellJacFactorsX.fill(static_cast<scalar_type>(1));
    bdCellJacFactorsY.fill(static_cast<scalar_type>(1));
    bdCellJacFactorsX[1] = static_cast<scalar_type>(-1);
    bdCellJacFactorsY[2] = static_cast<scalar_type>(-1);

    const auto & graphRows = m_meshObj.graphRowsOfCellsNearBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];

      FillStencilX(smPt, it, numDofPerCell);
      funcx(smPt, numDofPerCell, bdCellJacFactorsX);
      FillStencilY(smPt, it, numDofPerCell);
      funcy(smPt, numDofPerCell, bdCellJacFactorsY);
      addForcingContributionToVelocityAndJacobian(U, V, J, smPt);
    }
  }

  template<class U_t, class V_t>
  void velocityAndJacNearBDCellsImplDifferentScheme(const U_t & U,
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
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;

    // if here, then the velocity must be computed with Weno,
    /// while the jacobian must be computed with first order

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, stencil_container_type,
      U_t, MeshType, ghost_container_type>;

    // *****************************
    // *** functors for velocity ***
    // *****************************
    const auto stencilSizeForV = reconstructionTypeToStencilSize(m_recEn);
    stencil_container_type stencilValsForV(numDofPerCell*stencilSizeForV);

    stencil_filler_t FillStencilVeloX(reconstructionTypeToStencilSize(m_recEn),
				      U, m_meshObj, m_ghostLeft, m_ghostRight,
				      stencilValsForV, xAxis);
    stencil_filler_t FillStencilVeloY(reconstructionTypeToStencilSize(m_recEn),
				      U, m_meshObj, m_ghostBack, m_ghostFront,
				      stencilValsForV, yAxis);

    using velo_functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::implswe2d::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  scalar_type, flux_type>,
      V_t, scalar_type
      >;

    velo_functor_type funcVeloX(V, m_meshObj.dxInv(),
				/* end args for velo */
				m_fluxEn, normalX_, m_gravity, fluxL, fluxR,
				/* end args for flux */
				toReconstructionScheme(m_recEn), stencilValsForV,
				uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
				/* end args for reconstructor */
				);

    velo_functor_type funcVeloY(V, m_meshObj.dyInv(),
				/* end args for velo */
				m_fluxEn, normalY_, m_gravity, fluxB, fluxF,
				/* end args for flux */
				toReconstructionScheme(m_recEn), stencilValsForV,
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
				     U, m_meshObj, m_ghostLeft, m_ghostRight,
				     stencilValsForJ, xAxis);
    stencil_filler_t FillStencilJacY(stencilSizeForJ,
				     U, m_meshObj, m_ghostBack, m_ghostFront,
				     stencilValsForJ, yAxis);

    using jac_functor_type =
      pda::impl::ComputeDirectionalFluxBalanceFirstOrderJacobianOnBoundaryCell<
	pda::implswe2d::ComputeDirectionalFluxJacobians<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  scalar_type, flux_jac_type>,
      dimensionality, MeshType, jacobian_type
      >;

    jac_functor_type funcJacX(J, xAxis, m_meshObj,
			      /* end args for jac */
			      m_fluxEn, normalX_, m_gravity,
			      fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
			      /* end args for flux */
			      toReconstructionScheme(firstOrderRec), stencilValsForJ,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			      /* end args for reconstructor */
			      );

    jac_functor_type funcJacY(J, yAxis, m_meshObj,
			      /* end args for jac */
			      m_fluxEn, normalY_, m_gravity,
			      fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
			      /* end args for flux */
			      toReconstructionScheme(firstOrderRec), stencilValsForJ,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			      /* end args for reconstructor */
			      );

    std::array<scalar_type, numDofPerCell> bdCellJacFactorsX;
    std::array<scalar_type, numDofPerCell> bdCellJacFactorsY;
    bdCellJacFactorsX.fill(static_cast<scalar_type>(1));
    bdCellJacFactorsY.fill(static_cast<scalar_type>(1));
    bdCellJacFactorsX[1] = static_cast<scalar_type>(-1);
    bdCellJacFactorsY[2] = static_cast<scalar_type>(-1);

    // ************
    // loop
    // ************
    const auto & graphRows = m_meshObj.graphRowsOfCellsNearBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];
      FillStencilVeloX(smPt, it, numDofPerCell);
      funcVeloX(smPt, numDofPerCell);
      FillStencilJacX(smPt, it, numDofPerCell);
      funcJacX(smPt, numDofPerCell, bdCellJacFactorsX);

      FillStencilVeloY(smPt, it, numDofPerCell);
      funcVeloY(smPt, numDofPerCell);
      FillStencilJacY(smPt, it, numDofPerCell);
      funcJacY(smPt, numDofPerCell, bdCellJacFactorsY);

      addForcingContributionToVelocityAndJacobian(U, V, J, smPt);
    }
  }

  template<class U_t, class V_t>
  void velocityOnlyNearBdCellsImpl(const U_t & U,
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

    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    stencil_container_type stencilVals(numDofPerCell*stencilSize);

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, stencil_container_type, U_t, MeshType, ghost_container_type>;
    stencil_filler_t FillStencilX(reconstructionTypeToStencilSize(m_recEn),
				  U, m_meshObj, m_ghostLeft, m_ghostRight,
				  stencilVals, xAxis);

    stencil_filler_t FillStencilY(reconstructionTypeToStencilSize(m_recEn),
				  U, m_meshObj, m_ghostBack, m_ghostFront,
				  stencilVals, yAxis);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::implswe2d::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  scalar_type, flux_type>,
      V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.dxInv(),
		    /* end args for velo */
		    m_fluxEn, normalX_, m_gravity, fluxL, fluxR,
		    /* end args for flux */
		    toReconstructionScheme(m_recEn),stencilVals,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.dyInv(),
		    /* end args for velo */
		    m_fluxEn, normalY_, m_gravity, fluxB, fluxF,
		    /* end args for flux */
		    toReconstructionScheme(m_recEn),stencilVals,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		    /* end args for reconstructor */
		    );

    const auto & graphRows = m_meshObj.graphRowsOfCellsNearBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      FillStencilX(smPt, it, numDofPerCell);
      Fx(smPt, numDofPerCell);
      FillStencilY(smPt, it, numDofPerCell);
      Fy(smPt, numDofPerCell);
      addForcingContributionToVelocity(U, V, smPt);
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
	pda::implswe2d::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorForDiscreteFunction<
	    dimensionality, MeshType, U_t, edge_rec_type>,
	  scalar_type, flux_type>,
      V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.dxInv(),
		    /* end args for velo */
		    m_fluxEn, normalX_, m_gravity, fluxL, fluxR,
		    /* end args for flux */
		    xAxis, toReconstructionScheme(m_recEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.dyInv(),
		    /* end args for velo */
		    m_fluxEn, normalY_, m_gravity, fluxB, fluxF,
		    /* end args for flux */
		    yAxis, toReconstructionScheme(m_recEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		    /* end args for reconstructor */
		    );

    const auto & graphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      Fx(smPt, numDofPerCell);
      Fy(smPt, numDofPerCell);
      addForcingContributionToVelocity(U, V, smPt);
    }
  }

  template<class U_t, class V_t>
  void addForcingContributionToVelocity(const U_t & U,
					V_t & V,
					index_t smPt) const
  {
    const auto vIndex = smPt*numDofPerCell;
    const auto uIndex = m_meshObj.graph()(smPt, 0)*numDofPerCell;
    V(vIndex+1) -= m_coriolis*U(uIndex+2)/U(uIndex);
    V(vIndex+2) += m_coriolis*U(uIndex+1)/U(uIndex);
  }

  template<class U_t, class V_t>
  void addForcingContributionToVelocityAndJacobian(const U_t & U,
						   V_t & V,
						   jacobian_type & J,
						   index_t smPt) const
  {
    const auto vIndex = smPt*numDofPerCell;
    const index_t col_i = m_meshObj.graph()(smPt, 0)*numDofPerCell;
    V(vIndex+1) -= m_coriolis*U(col_i+2)/U(col_i);
    V(vIndex+2) += m_coriolis*U(col_i+1)/U(col_i);
    J.coeffRef(vIndex+1, col_i)   +=  m_coriolis*U(col_i+2)/(U(col_i)*U(col_i) );
    J.coeffRef(vIndex+1, col_i+2) += -m_coriolis/U(col_i);
    J.coeffRef(vIndex+2, col_i+1) +=  m_coriolis/U(col_i);
    J.coeffRef(vIndex+2, col_i)   +=  -m_coriolis*U(col_i+1)/(U(col_i)*U(col_i) );
  }

private:
  void allocateGhosts()
  {
    const auto stencilSize    = reconstructionTypeToStencilSize(m_recEn);
    const auto numGhostValues = numDofPerCell*((stencilSize-1)/2);

    const index_t s1 = m_meshObj.numCellsBd();
    ::pressiodemoapps::resize(m_ghostLeft, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostFront,s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight,s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostBack, s1, numGhostValues);
  }

protected:
  ::pressiodemoapps::Swe2d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  ::pressiodemoapps::InviscidFluxScheme m_fluxEn;
  const MeshType & m_meshObj;
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};
  mutable ghost_container_type m_ghostLeft;
  mutable ghost_container_type m_ghostFront;
  mutable ghost_container_type m_ghostRight;
  mutable ghost_container_type m_ghostBack;

  scalar_type m_gravity  = {};
  scalar_type m_coriolis = {};
  scalar_type m_initialPulseMagnitude = {};

  const std::array<scalar_type, 2> normalX_{1, 0};
  const std::array<scalar_type, 2> normalY_{0, 1};
};

template<class MeshType> constexpr int EigenApp<MeshType>::numDofPerCell;
template<class MeshType> constexpr int EigenApp<MeshType>::dimensionality;

}}//end namespace
#endif
