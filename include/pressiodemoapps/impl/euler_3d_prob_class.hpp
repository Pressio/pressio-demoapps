/*
//@HEADER
// ************************************************************************
//
// euler_3d_prob_class.hpp
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

#ifndef PRESSIODEMOAPPS_EULER3D_HPP_
#define PRESSIODEMOAPPS_EULER3D_HPP_

#include "euler_rankine_hugoniot.hpp"
#include "euler_rusanov_flux_values_function.hpp"
#include "euler_rusanov_flux_jacobian_function.hpp"
#include "euler_3d_initial_condition.hpp"
#include "euler_3d_ghost_filler_sedov.hpp"
#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "euler_flux_mixin.hpp"
#include "mixin_directional_flux_balance.hpp"
#include "mixin_directional_flux_balance_jacobian.hpp"

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#include <omp.h>
#endif

namespace pressiodemoapps{
namespace impleuler3d{

template<class MeshType>
class EigenApp
{

public:
  using index_t	      = typename MeshType::index_t;
  using scalar_type   = typename MeshType::scalar_t;
  using state_type    = Eigen::Matrix<scalar_type,-1,1>;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, index_t>;
  using mesh_connectivity_graph_type = typename MeshType::graph_t;

  static constexpr int dimensionality{3};
  static constexpr int numDofPerCell{5};

private:
  using ghost_container_type      = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using stencil_container_type    = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using flux_type	          = Eigen::Matrix<scalar_type, numDofPerCell,  1>;
  using edge_rec_type	          = Eigen::Matrix<scalar_type, numDofPerCell,  1>;
  using flux_jac_type             = Eigen::Matrix<scalar_type, numDofPerCell, numDofPerCell>;
  using reconstruction_gradient_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  EigenApp() = delete;

  EigenApp(const MeshType & meshObj,
	   ::pressiodemoapps::Euler3d probEnum,
	   ::pressiodemoapps::InviscidFluxReconstruction recEnum,
	   ::pressiodemoapps::InviscidFluxScheme fluxEnum)
    : m_meshObj(meshObj),
      m_probEn(probEnum),
      m_recEn(recEnum),
      m_fluxEn(fluxEnum)
  {
    m_numDofStencilMesh = m_meshObj.get().stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.get().sampleMeshSize() * numDofPerCell;
    allocateGhosts();
  }

  scalar_type gamma() const{
    return m_gamma;
  }

  state_type initialCondition() const
  {
    state_type IC(m_numDofStencilMesh);

    switch(m_probEn)
      {
      case ::pressiodemoapps::Euler3d::PeriodicSmooth:{
	euler3dsmoothInitialCondition(IC, m_meshObj.get(), m_gamma);
	return IC;
      }
      case ::pressiodemoapps::Euler3d::SedovSymmetry:{
	sedov3dInitialCondition(IC, m_meshObj.get(), m_gamma);
	return IC;
      }
      };

    return IC;
  }

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

    // reconstructions values
    edge_rec_type uMinusHalfNeg, uMinusHalfPos;
    edge_rec_type uPlusHalfNeg,  uPlusHalfPos;
    // fluxes
    flux_type fluxL, fluxR; //left, right x
    flux_type fluxB, fluxF; //back, front y
    flux_type fluxD, fluxU; //down, up z

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
			      fluxL, fluxR,
			      fluxB, fluxF,
			      fluxD, fluxU,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);
    }

    else{
      if (!m_meshObj.get().isFullyPeriodic()){
	velocityOnlyNearBdCellsImpl(U, currentTime, V,
				    fluxL, fluxR, fluxB, fluxF, fluxD, fluxU,
				    uMinusHalfNeg, uMinusHalfPos,
				    uPlusHalfNeg,  uPlusHalfPos);
      }

      velocityOnlyInnerCellsImpl(U, currentTime, V,
				 fluxL, fluxR, fluxB, fluxF, fluxD, fluxU,
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

    const auto zero = static_cast<scalar_type>(0);
    const auto & graph = m_meshObj.get().graph();
    const auto & targetGraphRows = m_meshObj.get().graphRowsOfCellsAwayFromBd();
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
	  (m_recEn == InviscidFluxReconstruction::FirstOrder) ? 6
	    : (m_recEn == InviscidFluxReconstruction::Weno3) ? 12 : 18;

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

    const auto zero = static_cast<scalar_type>(0);
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
	for (int i=1; i<=6; ++i){
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
    if (m_probEn == ::pressiodemoapps::Euler3d::SedovSymmetry)
    {
      using ghost_filler_t  = Ghost3dSedov<U_t, MeshType, ghost_container_type>;
      ghost_filler_t ghF(stencilSize, U, m_meshObj.get(),
			 m_ghostLeft, m_ghostRight,
			 m_ghostBack, m_ghostFront,
			 m_ghostDown, m_ghostUp);

      const auto & rowsBd = m_meshObj.get().graphRowsOfCellsNearBd();

      if (stencilSize==3){
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
	for (decltype(rowsBd.size()) it=0; it<rowsBd.size(); ++it){
	  ghF.template operator()<3>(rowsBd[it], it);
	}
      }
      else if (stencilSize==5){
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
	for (decltype(rowsBd.size()) it=0; it<rowsBd.size(); ++it){
	  ghF.template operator()<5>(rowsBd[it], it);
	}
      }else{
	throw std::runtime_error("missing impl");
      }
    }
  }

  template<class U_t, class V_t>
  void velocityAndJacobianImpl(const U_t & U,
			       const scalar_type currentTime,
			       V_t & V,
			       jacobian_type & J,
			       flux_type & fluxL,
			       flux_type & fluxR,
			       flux_type & fluxB,
			       flux_type & fluxF,
			       flux_type & fluxD,
			       flux_type & fluxU,
			       edge_rec_type & uMinusHalfNeg,
			       edge_rec_type & uMinusHalfPos,
			       edge_rec_type & uPlusHalfNeg,
			       edge_rec_type & uPlusHalfPos) const
  {

    // flux jacobians
    flux_jac_type fluxJacLNeg, fluxJacLPos;
    flux_jac_type fluxJacRNeg, fluxJacRPos;
    flux_jac_type fluxJacBNeg, fluxJacBPos;
    flux_jac_type fluxJacFNeg, fluxJacFPos;
    flux_jac_type fluxJacDNeg, fluxJacDPos;
    flux_jac_type fluxJacUNeg, fluxJacUPos;

    int nonZerosCountBeforeComputing = J.nonZeros();

    // near boundary I have be careful because
    // the jacobian can only be first order for now
    // only need to do near-BD cells if there are any
    if (!m_meshObj.get().isFullyPeriodic()){
      if (m_recEn == InviscidFluxReconstruction::FirstOrder){
	velocityAndJacNearBDCellsImplFirstOrder(U, currentTime, V, J,
						fluxL, fluxR,
						fluxB, fluxF,
						fluxD, fluxU,
						fluxJacLNeg, fluxJacLPos,
						fluxJacRNeg, fluxJacRPos,
						fluxJacBNeg, fluxJacBPos,
						fluxJacFNeg, fluxJacFPos,
						fluxJacDNeg, fluxJacDPos,
						fluxJacUNeg, fluxJacUPos,
						uMinusHalfNeg, uMinusHalfPos,
						uPlusHalfNeg,  uPlusHalfPos);
      }

      else{
	velocityAndJacNearBDCellsImplDifferentScheme(U, currentTime, V, J,
						     fluxL, fluxR,
						     fluxB, fluxF,
						     fluxD, fluxU,
						     fluxJacLNeg, fluxJacLPos,
						     fluxJacRNeg, fluxJacRPos,
						     fluxJacBNeg, fluxJacBPos,
						     fluxJacFNeg, fluxJacFPos,
						     fluxJacDNeg, fluxJacDPos,
						     fluxJacUNeg, fluxJacUPos,
						     uMinusHalfNeg, uMinusHalfPos,
						     uPlusHalfNeg,  uPlusHalfPos);
      }
    }

    velocityAndJacInnerCellsImpl(U, currentTime, V, J,
				 fluxL, fluxR,
				 fluxB, fluxF,
				 fluxD, fluxU,
				 fluxJacLNeg, fluxJacLPos,
				 fluxJacRNeg, fluxJacRPos,
				 fluxJacBNeg, fluxJacBPos,
				 fluxJacFNeg, fluxJacFPos,
				 fluxJacDNeg, fluxJacDPos,
				 fluxJacUNeg, fluxJacUPos,
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
				    flux_type & fluxR,
				    flux_type & fluxB,
				    flux_type & fluxF,
				    flux_type & fluxD,
				    flux_type & fluxU,
				    flux_jac_type & fluxJacLNeg,
				    flux_jac_type & fluxJacLPos,
				    flux_jac_type & fluxJacRNeg,
				    flux_jac_type & fluxJacRPos,
				    flux_jac_type & fluxJacBNeg,
				    flux_jac_type & fluxJacBPos,
				    flux_jac_type & fluxJacFNeg,
				    flux_jac_type & fluxJacFPos,
				    flux_jac_type & fluxJacDNeg,
				    flux_jac_type & fluxJacDPos,
				    flux_jac_type & fluxJacUNeg,
				    flux_jac_type & fluxJacUPos,
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
    constexpr int zAxis = 3;

    // allocate gradients of reconstructed states
    // the size depends on the scheme selected
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    reconstruction_gradient_t gradLNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradLPos(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRPos(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradBNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradBPos(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradFNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradFPos(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradDNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradDPos(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradUNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradUPos(numDofPerCell, stencilSize-1);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impl::ComputeDirectionalFluxBalanceJacobianOnInteriorCell<
	  pda::ee::impl::ComputeDirectionalFluxValuesAndJacobians<
	    pda::impl::ReconstructorForDiscreteFunction<
	      dimensionality, MeshType, U_t, edge_rec_type, reconstruction_gradient_t>,
	  numDofPerCell, scalar_type, flux_type, flux_jac_type>,
	dimensionality, MeshType, jacobian_type>,
      V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.get().dxInv(),
		    /* end args for velo */
		    J, xAxis, m_meshObj.get(),
		    /* end args for jac */
		    m_fluxEn, normalX_, m_gamma, fluxL, fluxR,
		    fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
		    /* end args for flux */
		    xAxis, toReconstructionScheme(m_recEn), U, m_meshObj.get(),
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		    gradLNeg, gradLPos, gradRNeg, gradRPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.get().dyInv(),
		    /* end args for velo */
		    J, yAxis, m_meshObj.get(),
		    /* end args for jac */
		    m_fluxEn, normalY_, m_gamma, fluxB, fluxF,
		    fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
		    /* end args for flux */
		    yAxis, toReconstructionScheme(m_recEn), U, m_meshObj.get(),
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		    gradBNeg, gradBPos, gradFNeg, gradFPos
		    /* end args for reconstructor */
		    );

    functor_type Fz(V, m_meshObj.get().dzInv(),
		    /* end args for velo */
		    J, zAxis, m_meshObj.get(),
		    /* end args for jac */
		    m_fluxEn, normalZ_, m_gamma, fluxD, fluxU,
		    fluxJacDNeg, fluxJacDPos, fluxJacUNeg, fluxJacUPos,
		    /* end args for flux */
		    zAxis, toReconstructionScheme(m_recEn), U, m_meshObj.get(),
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		    gradDNeg, gradDPos, gradUNeg, gradUPos
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
      Fz(smPt, numDofPerCell);
    }
  }

  template<class U_t, class V_t>
  void velocityAndJacNearBDCellsImplFirstOrder(const U_t & U,
					       const scalar_type /*currentTime*/,
					       V_t & V,
					       jacobian_type & J,
					       flux_type & fluxL,
					       flux_type & fluxR,
					       flux_type & fluxB,
					       flux_type & fluxF,
					       flux_type & fluxD,
					       flux_type & fluxU,
					       flux_jac_type & fluxJacLNeg,
					       flux_jac_type & fluxJacLPos,
					       flux_jac_type & fluxJacRNeg,
					       flux_jac_type & fluxJacRPos,
					       flux_jac_type & fluxJacBNeg,
					       flux_jac_type & fluxJacBPos,
					       flux_jac_type & fluxJacFNeg,
					       flux_jac_type & fluxJacFPos,
					       flux_jac_type & fluxJacDNeg,
					       flux_jac_type & fluxJacDPos,
					       flux_jac_type & fluxJacUNeg,
					       flux_jac_type & fluxJacUPos,
					       edge_rec_type & uMinusHalfNeg,
					       edge_rec_type & uMinusHalfPos,
					       edge_rec_type & uPlusHalfNeg,
					       edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    constexpr int zAxis = 3;
    assert(m_recEn == InviscidFluxReconstruction::FirstOrder);

    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    stencil_container_type stencilVals(numDofPerCell*stencilSize);

    // if here, then the scheme for velocity matches
    // the one for Jacobian so we can use same functors

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, stencil_container_type, U_t, MeshType, ghost_container_type>;
    stencil_filler_t FillStencilX(reconstructionTypeToStencilSize(m_recEn),
				   U, m_meshObj.get(), m_ghostLeft, m_ghostRight,
				  stencilVals, xAxis);

    stencil_filler_t FillStencilY(reconstructionTypeToStencilSize(m_recEn),
				   U, m_meshObj.get(), m_ghostBack, m_ghostFront,
				   stencilVals, yAxis);

    stencil_filler_t FillStencilZ(reconstructionTypeToStencilSize(m_recEn),
				   U, m_meshObj.get(), m_ghostDown, m_ghostUp,
				   stencilVals, zAxis);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impl::ComputeDirectionalFluxBalanceFirstOrderJacobianOnBoundaryCell<
	  pda::ee::impl::ComputeDirectionalFluxValuesAndJacobians<
	    pda::impl::ReconstructorFromStencil<
	      edge_rec_type, stencil_container_type>,
	    numDofPerCell, scalar_type, flux_type, flux_jac_type>,
	  dimensionality, MeshType, jacobian_type>,
      V_t, scalar_type
      >;

    functor_type funcx(V, m_meshObj.get().dxInv(),
		       /* end args for velo */
		       J, xAxis, m_meshObj.get(),
		       /* end args for jac */
		       m_fluxEn, normalX_, m_gamma, fluxL, fluxR,
		       fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
		       /* end args for flux */
		       toReconstructionScheme(m_recEn), stencilVals,
		       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		       /* end args for reconstructor */
		       );

    functor_type funcy(V, m_meshObj.get().dyInv(),
		       /* end args for velo */
		       J, yAxis, m_meshObj.get(),
		       /* end args for jac */
		       m_fluxEn, normalY_, m_gamma, fluxB, fluxF,
		       fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
		       /* end args for flux */
		       toReconstructionScheme(m_recEn), stencilVals,
		       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		       /* end args for reconstructor */
		       );

    functor_type funcz(V, m_meshObj.get().dzInv(),
		       /* end args for velo */
		       J, zAxis, m_meshObj.get(),
		       /* end args for jac */
		       m_fluxEn, normalZ_, m_gamma, fluxD, fluxU,
		       fluxJacDNeg, fluxJacDPos, fluxJacUNeg, fluxJacUPos,
		       /* end args for flux */
		       toReconstructionScheme(m_recEn), stencilVals,
		       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		       /* end args for reconstructor */
		       );

    std::array<scalar_type, numDofPerCell> bcCellJacFactorsDefault;
    bcCellJacFactorsDefault.fill(static_cast<scalar_type>(1));

    std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveX;
    std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveY;
    std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveZ;

    bcCellJacFactorsReflectiveX.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsReflectiveY.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsReflectiveZ.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsReflectiveX[1] = static_cast<scalar_type>(-1);
    bcCellJacFactorsReflectiveY[2] = static_cast<scalar_type>(-1);
    bcCellJacFactorsReflectiveZ[3] = static_cast<scalar_type>(-1);

    const auto & graphRows = m_meshObj.get().graphRowsOfCellsNearBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];

      FillStencilX(smPt, it, numDofPerCell);
      auto bcTypeX = findCellBdType(smPt, xAxis);
      const auto & factorsX = (bcTypeX == 1)
	? bcCellJacFactorsReflectiveX : bcCellJacFactorsDefault;
      funcx(smPt, numDofPerCell, factorsX);

      FillStencilY(smPt, it, numDofPerCell);
      auto bcTypeY = findCellBdType(smPt, yAxis);
      const auto & factorsY = (bcTypeY == 1)
	? bcCellJacFactorsReflectiveY : bcCellJacFactorsDefault;
      funcy(smPt, numDofPerCell, factorsY);

      FillStencilZ(smPt, it, numDofPerCell);
      auto bcTypeZ = findCellBdType(smPt, zAxis);
      const auto & factorsZ = (bcTypeZ == 1)
	? bcCellJacFactorsReflectiveZ : bcCellJacFactorsDefault;
      funcz(smPt, numDofPerCell, factorsZ);
    }
  }

  template<class U_t, class V_t>
  void velocityAndJacNearBDCellsImplDifferentScheme(const U_t & U,
						    const scalar_type /*currentTime*/,
						    V_t & V,
						    jacobian_type & J,
						    flux_type & fluxL,
						    flux_type & fluxR,
						    flux_type & fluxB,
						    flux_type & fluxF,
						    flux_type & fluxD,
						    flux_type & fluxU,
						    flux_jac_type & fluxJacLNeg,
						    flux_jac_type & fluxJacLPos,
						    flux_jac_type & fluxJacRNeg,
						    flux_jac_type & fluxJacRPos,
						    flux_jac_type & fluxJacBNeg,
						    flux_jac_type & fluxJacBPos,
						    flux_jac_type & fluxJacFNeg,
						    flux_jac_type & fluxJacFPos,
						    flux_jac_type & fluxJacDNeg,
						    flux_jac_type & fluxJacDPos,
						    flux_jac_type & fluxJacUNeg,
						    flux_jac_type & fluxJacUPos,
						    edge_rec_type & uMinusHalfNeg,
						    edge_rec_type & uMinusHalfPos,
						    edge_rec_type & uPlusHalfNeg,
						    edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    constexpr int zAxis = 3;

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
				      U, m_meshObj.get(), m_ghostLeft, m_ghostRight,
				      stencilValsForV, xAxis);
    stencil_filler_t FillStencilVeloY(reconstructionTypeToStencilSize(m_recEn),
				      U, m_meshObj.get(), m_ghostBack, m_ghostFront,
				      stencilValsForV, yAxis);
    stencil_filler_t FillStencilVeloZ(reconstructionTypeToStencilSize(m_recEn),
				      U, m_meshObj.get(), m_ghostDown, m_ghostUp,
				      stencilValsForV, zAxis);

    using velo_functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::ee::impl::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  numDofPerCell, scalar_type, flux_type>,
      V_t, scalar_type
      >;

    velo_functor_type funcVeloX(V, m_meshObj.get().dxInv(),
				/* end args for velo */
				m_fluxEn, normalX_, m_gamma, fluxL, fluxR,
				/* end args for flux */
				toReconstructionScheme(m_recEn), stencilValsForV,
				uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
				/* end args for reconstructor */
				);

    velo_functor_type funcVeloY(V, m_meshObj.get().dyInv(),
				/* end args for velo */
				m_fluxEn, normalY_, m_gamma, fluxB, fluxF,
				/* end args for flux */
				toReconstructionScheme(m_recEn), stencilValsForV,
				uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
				/* end args for reconstructor */
				);

    velo_functor_type funcVeloZ(V, m_meshObj.get().dzInv(),
				/* end args for velo */
				m_fluxEn, normalZ_, m_gamma, fluxD, fluxU,
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
				     U, m_meshObj.get(), m_ghostLeft, m_ghostRight,
				     stencilValsForJ, xAxis);
    stencil_filler_t FillStencilJacY(stencilSizeForJ,
				     U, m_meshObj.get(), m_ghostBack, m_ghostFront,
				     stencilValsForJ, yAxis);
    stencil_filler_t FillStencilJacZ(stencilSizeForJ,
				     U, m_meshObj.get(), m_ghostDown, m_ghostUp,
				     stencilValsForJ, zAxis);

    using jac_functor_type =
      pda::impl::ComputeDirectionalFluxBalanceFirstOrderJacobianOnBoundaryCell<
	pda::ee::impl::ComputeDirectionalFluxJacobians<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  numDofPerCell, scalar_type, flux_jac_type>,
      dimensionality, MeshType, jacobian_type
      >;

    jac_functor_type funcJacX(J, xAxis, m_meshObj.get(),
			      /* end args for jac */
			      m_fluxEn, normalX_, m_gamma,
			      fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
			      /* end args for flux */
			      toReconstructionScheme(firstOrderRec), stencilValsForJ,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			      /* end args for reconstructor */
			      );

    jac_functor_type funcJacY(J, yAxis, m_meshObj.get(),
			      /* end args for jac */
			      m_fluxEn, normalY_, m_gamma,
			      fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
			      /* end args for flux */
			      toReconstructionScheme(firstOrderRec), stencilValsForJ,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			      /* end args for reconstructor */
			      );

    jac_functor_type funcJacZ(J, zAxis, m_meshObj.get(),
			      /* end args for jac */
			      m_fluxEn, normalZ_, m_gamma,
			      fluxJacDNeg, fluxJacDPos, fluxJacUNeg, fluxJacUPos,
			      /* end args for flux */
			      toReconstructionScheme(firstOrderRec), stencilValsForJ,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			      /* end args for reconstructor */
			      );

    std::array<scalar_type, numDofPerCell> bcCellJacFactorsDefault;
    bcCellJacFactorsDefault.fill(static_cast<scalar_type>(1));

    std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveX;
    std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveY;
    std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveZ;
    bcCellJacFactorsReflectiveX.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsReflectiveY.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsReflectiveZ.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsReflectiveX[1] = static_cast<scalar_type>(-1);
    bcCellJacFactorsReflectiveY[2] = static_cast<scalar_type>(-1);
    bcCellJacFactorsReflectiveZ[3] = static_cast<scalar_type>(-1);

    // ************
    // loop
    // ************
    const auto & graphRows = m_meshObj.get().graphRowsOfCellsNearBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      FillStencilVeloX(smPt, it, numDofPerCell);
      funcVeloX(smPt, numDofPerCell);
      FillStencilJacX(smPt, it, numDofPerCell);
      auto bcTypeX = findCellBdType(smPt, xAxis);
      const auto & factorsX = (bcTypeX == 1) ? bcCellJacFactorsReflectiveX : bcCellJacFactorsDefault;
      funcJacX(smPt, numDofPerCell, factorsX);

      FillStencilVeloY(smPt, it, numDofPerCell);
      funcVeloY(smPt, numDofPerCell);
      FillStencilJacY(smPt, it, numDofPerCell);
      auto bcTypeY = findCellBdType(smPt, yAxis);
      const auto & factorsY = (bcTypeY == 1) ? bcCellJacFactorsReflectiveY : bcCellJacFactorsDefault;
      funcJacY(smPt, numDofPerCell, factorsY);

      FillStencilVeloZ(smPt, it, numDofPerCell);
      funcVeloZ(smPt, numDofPerCell);
      FillStencilJacZ(smPt, it, numDofPerCell);
      auto bcTypeZ = findCellBdType(smPt, zAxis);
      const auto & factorsZ = (bcTypeZ == 1) ? bcCellJacFactorsReflectiveZ : bcCellJacFactorsDefault;
      funcJacZ(smPt, numDofPerCell, factorsZ);

    }
  }

  template<class U_t, class V_t>
  void velocityOnlyInnerCellsImpl(const U_t & U,
				  const scalar_type /*currentTime*/,
				  V_t & V,
				  flux_type & fluxL,
				  flux_type & fluxR,
				  flux_type & fluxB,
				  flux_type & fluxF,
				  flux_type & fluxD,
				  flux_type & fluxU,
				  edge_rec_type & uMinusHalfNeg,
				  edge_rec_type & uMinusHalfPos,
				  edge_rec_type & uPlusHalfNeg,
				  edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    constexpr int zAxis = 3;

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::ee::impl::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorForDiscreteFunction<
	    dimensionality, MeshType, U_t, edge_rec_type>,
	  numDofPerCell, scalar_type, flux_type>,
      V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.get().dxInv(),
		    /* end args for velo */
		    m_fluxEn, normalX_, m_gamma, fluxL, fluxR,
		    /* end args for flux */
		    xAxis, toReconstructionScheme(m_recEn), U, m_meshObj.get(),
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.get().dyInv(),
		    /* end args for velo */
		    m_fluxEn, normalY_, m_gamma, fluxB, fluxF,
		    /* end args for flux */
		    yAxis, toReconstructionScheme(m_recEn), U, m_meshObj.get(),
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fz(V, m_meshObj.get().dzInv(),
		    /* end args for velo */
		    m_fluxEn, normalZ_, m_gamma, fluxD, fluxU,
		    /* end args for flux */
		    zAxis, toReconstructionScheme(m_recEn), U, m_meshObj.get(),
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
      Fz(smPt, numDofPerCell);
    }
  }

  template<class U_t, class V_t>
  void velocityOnlyNearBdCellsImpl(const U_t & U,
				   const scalar_type /*currentTime*/,
				   V_t & V,
				   flux_type & fluxL,
				   flux_type & fluxR,
				   flux_type & fluxB,
				   flux_type & fluxF,
				   flux_type & fluxD,
				   flux_type & fluxU,
				   edge_rec_type & uMinusHalfNeg,
				   edge_rec_type & uMinusHalfPos,
				   edge_rec_type & uPlusHalfNeg,
				   edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    constexpr int zAxis = 3;

    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    stencil_container_type stencilVals(numDofPerCell*stencilSize);

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, stencil_container_type, U_t, MeshType, ghost_container_type>;
    stencil_filler_t FillStencilX(reconstructionTypeToStencilSize(m_recEn),
				  U, m_meshObj.get(), m_ghostLeft, m_ghostRight,
				  stencilVals, xAxis);

    stencil_filler_t FillStencilY(reconstructionTypeToStencilSize(m_recEn),
				  U, m_meshObj.get(), m_ghostBack, m_ghostFront,
				  stencilVals, yAxis);

    stencil_filler_t FillStencilZ(reconstructionTypeToStencilSize(m_recEn),
				  U, m_meshObj.get(), m_ghostDown, m_ghostUp,
				  stencilVals, zAxis);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::ee::impl::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  numDofPerCell, scalar_type, flux_type>,
      V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.get().dxInv(),
		    /* end args for velo */
		    m_fluxEn, normalX_, m_gamma, fluxL, fluxR,
		    /* end args for flux */
		    toReconstructionScheme(m_recEn), stencilVals,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.get().dyInv(),
		    /* end args for velo */
		    m_fluxEn, normalY_, m_gamma, fluxB, fluxF,
		    /* end args for flux */
		    toReconstructionScheme(m_recEn), stencilVals,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fz(V, m_meshObj.get().dzInv(),
		    /* end args for velo */
		    m_fluxEn, normalZ_, m_gamma, fluxD, fluxU,
		    /* end args for flux */
		    toReconstructionScheme(m_recEn), stencilVals,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		    /* end args for reconstructor */
		    );

    const auto & graphRows = m_meshObj.get().graphRowsOfCellsNearBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      FillStencilX(smPt, it, numDofPerCell);
      Fx(smPt, numDofPerCell);
      FillStencilY(smPt, it, numDofPerCell);
      Fy(smPt, numDofPerCell);
      FillStencilZ(smPt, it, numDofPerCell);
      Fz(smPt, numDofPerCell);
    }
  }

  int findCellBdType(index_t graphRow, int axis) const
  {
    // 0: Neumann
    // 1: Reflective
    // 2: Dirichlet
    constexpr int neumann = 0;
    constexpr int reflective = 1;
    // constexpr int dirichlet = 2;

    if (m_probEn == ::pressiodemoapps::Euler3d::SedovSymmetry)
    {
      if (axis == 1 && m_meshObj.get().hasBdLeft3d(graphRow)){
	return reflective;
      }

      if (axis == 2 && m_meshObj.get().hasBdBack3d(graphRow)){
	return reflective;
      }

      if (axis == 3 && m_meshObj.get().hasBdBottom3d(graphRow)){
	return reflective;
      }

      return neumann;
    }

    return 0;
  }

  void allocateGhosts()
  {
    const auto stencilSize    = reconstructionTypeToStencilSize(m_recEn);
    const auto numGhostValues = numDofPerCell*((stencilSize-1)/2);

    const index_t s1 = m_meshObj.get().numCellsNearBd();
    ::pressiodemoapps::resize(m_ghostLeft,  s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostFront, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostBack,  s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostDown, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostUp,  s1, numGhostValues);
  }

protected:
  scalar_type m_gamma = static_cast<scalar_type>(1.4);

  std::reference_wrapper<const MeshType> m_meshObj;
  ::pressiodemoapps::Euler3d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  ::pressiodemoapps::InviscidFluxScheme m_fluxEn;

  // note that dof refers to the degress of freedom,
  // which is NOT same as grid points.
  // SampleMesh_ identifies the velocity/residual locations
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  mutable ghost_container_type m_ghostLeft;
  mutable ghost_container_type m_ghostRight;
  mutable ghost_container_type m_ghostBack;
  mutable ghost_container_type m_ghostFront;
  mutable ghost_container_type m_ghostDown;
  mutable ghost_container_type m_ghostUp;

  std::array<scalar_type, dimensionality> normalX_{1, 0, 0};
  std::array<scalar_type, dimensionality> normalY_{0, 1, 0};
  std::array<scalar_type, dimensionality> normalZ_{0, 0, 1};
};

template<class MeshType> constexpr int EigenApp<MeshType>::numDofPerCell;
template<class MeshType> constexpr int EigenApp<MeshType>::dimensionality;

}}//end namespace pressiodemoapps::ee::impl
#endif
