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
#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "advection_diffusion_2d_flux_mixin.hpp"
#include "mixin_directional_flux_balance.hpp"
#include "mixin_directional_flux_balance_jacobian.hpp"
#include "Eigen/Sparse"

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
  // constructor for Burgers2d periodic
  //
  EigenApp(TagBurgersPeriodic /*tag*/,
	   const MeshType & meshObj,
	   ::pressiodemoapps::InviscidFluxReconstruction inviscidFluxRecEn,
	   ::pressiodemoapps::InviscidFluxScheme invFluxSchemeEn,
	   ::pressiodemoapps::ViscousFluxReconstruction visFluxRecEn,
	   scalar_type icPulseMagnitude,
	   scalar_type icSpread,
	   scalar_type diffusionCoeff,
	   scalar_type x0,
	   scalar_type y0)
    : m_numDofPerCell(2),
      m_probEn(::pressiodemoapps::AdvectionDiffusion2d::BurgersPeriodic),
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
    m_numDofStencilMesh = m_meshObj.stencilMeshSize() * m_numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * m_numDofPerCell;
    // don't need to allocate ghosts because it is periodic
  }

  state_type initialCondition() const
  {
    state_type initialState(m_numDofStencilMesh);

    if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersPeriodic)
      {
	burgers2d_gaussian(initialState, m_meshObj,
			   m_burgers2d_icPulse,
			   m_burgers2d_icSpread,
			   m_burgers2d_x0,
			   m_burgers2d_y0);
      }

    return initialState;
  }

protected:
  void initializeJacobian(jacobian_type & J)
  {
    J.resize(m_numDofSampleMesh, m_numDofStencilMesh);

    using Tr = Eigen::Triplet<scalar_type>;
    std::vector<Tr> trList;

    // only need this since this is for now periodic
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

      if (J){
	velocityAndJacobianPeriodicImpl(U, currentTime, V, *J,
					fluxL, fluxF, fluxR, fluxB,
					uMinusHalfNeg, uMinusHalfPos,
					uPlusHalfNeg,  uPlusHalfPos);
      }
      else{

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

  template<class U_t, class V_t>
  void velocityAndJacobianPeriodicImpl(const U_t & U,
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

    // this is for periodic so we don't have boundary cells,
    // all cells count as "inner cells"
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
	  pda::impladvdiff2d::ComputeDirectionalFluxValuesAndJacobians<
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
		    m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR,
		    fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
		    /* end args for flux */
		    xAxis, toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		    gradLNeg, gradLPos, gradRNeg, gradRPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.dyInv(),
		    /* end args for velo */
		    J, yAxis, m_meshObj,
		    /* end args for jac */
		    m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF,
		    fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
		    /* end args for flux */
		    yAxis, toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj,
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
      Fx(smPt, m_numDofPerCell);
      Fy(smPt, m_numDofPerCell);

      if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersPeriodic){
	addBurgersDiffusionAndSourceToVelocityAndJacobianInnerCells(U, V, J, smPt);
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

    functor_type Fx(V, m_meshObj.dxInv(),
		    /* end args for velo */
		    m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR,
		    /* end args for flux */
		    xAxis, toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.dyInv(),
		    /* end args for velo */
		    m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF,
		    /* end args for flux */
		    yAxis, toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		    /* end args for reconstructor */
		    );

    const auto & graphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      Fx(smPt, m_numDofPerCell);
      Fy(smPt, m_numDofPerCell);

      if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersPeriodic){
	addBurgersDiffusionToVelocityInnerCells(U, V, smPt);
      }
    }
  }

  template<class U_t, class V_t>
  void addBurgersDiffusionToVelocityInnerCells(const U_t & U,
					       V_t & V,
					       index_t smPt) const
  {
    constexpr auto two  = static_cast<scalar_type>(2);
    const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto dyInvSq  = m_meshObj.dyInv()*m_meshObj.dyInv();
    const auto diffDxInvSq  = m_burgers2d_diffusion*dxInvSq;
    const auto diffDyInvSq  = m_burgers2d_diffusion*dyInvSq;

    const auto & graph     = m_meshObj.graph();
    const auto vIndex      = smPt*m_numDofPerCell;
    const auto uIndex      = graph(smPt, 0)*m_numDofPerCell;
    const auto uIndexLeft  = graph(smPt, 1)*m_numDofPerCell;
    const auto uIndexFront = graph(smPt, 2)*m_numDofPerCell;
    const auto uIndexRight = graph(smPt, 3)*m_numDofPerCell;
    const auto uIndexBack  = graph(smPt, 4)*m_numDofPerCell;
    V(vIndex) += diffDxInvSq*( U(uIndexRight) - two*U(uIndex) + U(uIndexLeft) );
    V(vIndex) += diffDyInvSq*( U(uIndexFront) - two*U(uIndex) + U(uIndexBack) );

    V(vIndex+1) += diffDxInvSq*( U(uIndexRight+1) - two*U(uIndex+1) + U(uIndexLeft+1) );
    V(vIndex+1) += diffDyInvSq*( U(uIndexFront+1) - two*U(uIndex+1) + U(uIndexBack+1) );
  }

  template<class U_t, class V_t>
  void addBurgersDiffusionAndSourceToVelocityAndJacobianInnerCells(const U_t & U,
								   V_t & V,
								   jacobian_type & J,
								   index_t smPt) const
  {
    constexpr auto two  = static_cast<scalar_type>(2);
    const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto dyInvSq  = m_meshObj.dyInv()*m_meshObj.dyInv();
    const auto diffDxInvSq  = m_burgers2d_diffusion*dxInvSq;
    const auto diffDyInvSq  = m_burgers2d_diffusion*dyInvSq;

    const auto & graph     = m_meshObj.graph();
    const auto vIndex      = smPt*m_numDofPerCell;
    const auto uIndex      = graph(smPt, 0)*m_numDofPerCell;
    const auto uIndexLeft  = graph(smPt, 1)*m_numDofPerCell;
    const auto uIndexFront = graph(smPt, 2)*m_numDofPerCell;
    const auto uIndexRight = graph(smPt, 3)*m_numDofPerCell;
    const auto uIndexBack  = graph(smPt, 4)*m_numDofPerCell;
    V(vIndex) += diffDxInvSq*( U(uIndexRight) - two*U(uIndex) + U(uIndexLeft) );
    V(vIndex) += diffDyInvSq*( U(uIndexFront) - two*U(uIndex) + U(uIndexBack) );

    V(vIndex+1) += diffDxInvSq*( U(uIndexRight+1) - two*U(uIndex+1) + U(uIndexLeft+1) );
    V(vIndex+1) += diffDyInvSq*( U(uIndexFront+1) - two*U(uIndex+1) + U(uIndexBack+1) );

    J.coeffRef(vIndex, uIndex) += -two*diffDxInvSq - two*diffDyInvSq;
    J.coeffRef(vIndex, uIndexLeft)  += diffDxInvSq;
    J.coeffRef(vIndex, uIndexFront) += diffDyInvSq;
    J.coeffRef(vIndex, uIndexRight) += diffDxInvSq;
    J.coeffRef(vIndex, uIndexBack)  += diffDyInvSq;
    J.coeffRef(vIndex+1, uIndex+1) += -two*diffDxInvSq - two*diffDyInvSq;
    J.coeffRef(vIndex+1, uIndexLeft+1)  += diffDxInvSq;
    J.coeffRef(vIndex+1, uIndexFront+1) += diffDyInvSq;
    J.coeffRef(vIndex+1, uIndexRight+1) += diffDxInvSq;
    J.coeffRef(vIndex+1, uIndexBack+1)  += diffDyInvSq;
  }

protected:
  // common to all problems
  int m_numDofPerCell = {};
  ::pressiodemoapps::AdvectionDiffusion2d m_probEn;
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
  scalar_type m_burgers2d_icPulse = {};
  scalar_type m_burgers2d_icSpread = {};
  scalar_type m_burgers2d_diffusion = {};
  scalar_type m_burgers2d_x0 = {};
  scalar_type m_burgers2d_y0 = {};
};

template<class MeshType> constexpr int EigenApp<MeshType>::dimensionality;

}}//end namespace
#endif
