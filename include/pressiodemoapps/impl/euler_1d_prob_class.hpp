/*
//@HEADER
// ************************************************************************
//
// euler_1d_prob_class.hpp
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

#ifndef PRESSIODEMOAPPS_EULER1D_APP_HPP_
#define PRESSIODEMOAPPS_EULER1D_APP_HPP_

#include "euler_rusanov_flux_values_function.hpp"
#include "euler_rusanov_flux_jacobian_function.hpp"
#include "euler_rankine_hugoniot.hpp"
#include "euler_1d_initial_condition.hpp"
#include "euler_1d_ghost_filler.hpp"
#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "euler_flux_mixin.hpp"
#include "mixin_directional_flux_balance.hpp"
#include "mixin_directional_flux_balance_jacobian.hpp"
#include "Eigen/Sparse"

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#include <omp.h>
#endif

namespace pressiodemoapps{ namespace impleuler1d{

template<class MeshType>
class EigenApp
{
public:
  // required public aliases
  using index_t	      = typename MeshType::index_t;
  using scalar_type   = typename MeshType::scalar_t;
  using state_type    = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, index_t>;
  using mesh_connectivity_graph_type = typename MeshType::graph_t;

private:
  static constexpr int dimensionality{1};
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
  EigenApp() = delete;

  EigenApp(const MeshType & meshObj,
	   ::pressiodemoapps::Euler1d probEnum,
	   ::pressiodemoapps::InviscidFluxReconstruction recEnum,
	   ::pressiodemoapps::InviscidFluxScheme fluxEnum)
    : m_numDofPerCell(3),
      m_probEn(probEnum),
      m_inviscidFluxRecEn(recEnum),
      m_inviscidFluxSchemeEn(fluxEnum),
      m_gamma(static_cast<scalar_type>(1.4)),
      m_meshObj(meshObj)
  {

    m_numDofStencilMesh = m_meshObj.get().stencilMeshSize() * m_numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.get().sampleMeshSize() * m_numDofPerCell;
    allocateGhostValues();
  }

  scalar_type gamma() const{ return m_gamma; }

  state_type initialCondition() const
  {
    state_type initialState(m_numDofStencilMesh);

    if (m_probEn == ::pressiodemoapps::Euler1d::PeriodicSmooth){
      euler1dsineInitialCondition(initialState, m_meshObj.get(), m_gamma);
    }
    else if (m_probEn == ::pressiodemoapps::Euler1d::Sod){
      sod1dInitialCondition(initialState, m_meshObj.get(), m_gamma);
    }
    else if (m_probEn == ::pressiodemoapps::Euler1d::Lax){
      lax1dInitialCondition(initialState, m_meshObj.get(), m_gamma);
    }
    else if (m_probEn == ::pressiodemoapps::Euler1d::ShuOsher){
      shuOsherInitialCondition(initialState, m_meshObj.get(), m_gamma);
    }
    else{
      //nothing
    }
    return initialState;
  }

protected:
  int numDofPerCellImpl() const {
    return m_numDofPerCell;
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

    // this is 1d, there is no point in fillin ghosts in parallel
    // so do this before starting parallel region
    fillGhostsIfNeeded(U);

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp parallel
    {
#endif

      // reconstructions values
      edge_rec_type uMinusHalfNeg(m_numDofPerCell);
      edge_rec_type uMinusHalfPos(m_numDofPerCell);
      edge_rec_type uPlusHalfNeg(m_numDofPerCell);
      edge_rec_type uPlusHalfPos(m_numDofPerCell);
      // fluxes
      flux_type fluxL(m_numDofPerCell);
      flux_type fluxR(m_numDofPerCell);

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

	velocityAndJacobianImpl(U, currentTime, V, *J,
				fluxL, fluxR,
				uMinusHalfNeg, uMinusHalfPos,
				uPlusHalfNeg,  uPlusHalfPos);
      }
      else
      {
	if (!m_meshObj.get().isFullyPeriodic()){
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp single
#endif
	  {
	    velocityOnlyNearBdCellsImpl(U, currentTime, V, fluxL, fluxR,
					uMinusHalfNeg, uMinusHalfPos,
					uPlusHalfNeg,  uPlusHalfPos);
	  }
	}

	velocityOnlyInnerCellsImpl(U, currentTime, V, fluxL, fluxR,
				   uMinusHalfNeg, uMinusHalfPos,
				   uPlusHalfNeg,  uPlusHalfPos);
      }

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
 }
#endif

  }

private:
  template<class Tr>
  void initializeJacobianForInnerCells(std::vector<Tr> & trList)
  {
    // for inner cells, the Jacobian is of the same scheme
    // wanted by the user, no special treatment is needed

    const scalar_type zero = 0;
    const auto & graph = m_meshObj.get().graph();
    const auto & targetGraphRows = m_meshObj.get().graphRowsOfCellsAwayFromBd();
    for (std::size_t it=0; it<targetGraphRows.size(); ++it)
      {
	const auto smPt = targetGraphRows[it];
	const auto jacRowOfCurrCellRho = smPt*m_numDofPerCell;

	// entries wrt current cell's dofs
	const auto jacColOfCurrCellRho = graph(smPt, 0)*m_numDofPerCell;
	for (int k=0; k<m_numDofPerCell; ++k){
	  for (int j=0; j<m_numDofPerCell; ++j){
	    trList.push_back( Tr(jacRowOfCurrCellRho+k, jacColOfCurrCellRho+j, zero) );
	  }
	}

	// wrt neighbors: this depends on the advection scheme
	const auto numNeighbors =
	  (m_inviscidFluxRecEn == InviscidFluxReconstruction::FirstOrder) ? 2
	    : (m_inviscidFluxRecEn == InviscidFluxReconstruction::Weno3) ? 4 : 6;

	for (int i=1; i<=numNeighbors; ++i){
	  const auto ci = graph(smPt, i)*m_numDofPerCell;
	  for (int k=0; k<m_numDofPerCell; ++k){
	    for (int j=0; j<m_numDofPerCell; ++j){
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
    const auto & graph = m_meshObj.get().graph();
    const auto & targetGraphRows = m_meshObj.get().graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<targetGraphRows.size(); ++it)
      {
	const auto smPt = targetGraphRows[it];
	const auto jacRowOfCurrCellRho = smPt*m_numDofPerCell;
	const auto jacColOfCurrCellRho = graph(smPt, 0)*m_numDofPerCell;

	// wrt current cell's dofs
	for (int k=0; k<m_numDofPerCell; ++k){
	  for (int j=0; j<m_numDofPerCell; ++j){
	    trList.push_back( Tr(jacRowOfCurrCellRho+k,
				 jacColOfCurrCellRho+j,
				 zero) );
	  }
	}

	// wrt neighbors
	// for near-bd, we only do first-order Jacobian for now
	const auto L0 = graph(smPt, 1);
	if (L0 != -1){
	  const auto ciL0 = L0*m_numDofPerCell;
	  for (int k=0; k<m_numDofPerCell; ++k){
	    for (int j=0; j<m_numDofPerCell; ++j){
	      trList.push_back( Tr(jacRowOfCurrCellRho+k, ciL0+j, zero) );
	    }
	  }
	}

	const auto R0 = graph(smPt, 2);
	if (R0 != -1){
	  const auto ciR0 = R0*m_numDofPerCell;
	  for (int k=0; k<m_numDofPerCell; ++k){
	    for (int j=0; j<m_numDofPerCell; ++j){
	      trList.push_back( Tr(jacRowOfCurrCellRho+k, ciR0+j, zero) );
	    }
	  }
	}
      }
  }

  template<class U_t>
  void fillGhostsIfNeeded(const U_t & U) const
  {
    namespace pda = ::pressiodemoapps;

    const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    // only need ghosts for specific problems
    if (m_probEn == pda::Euler1d::Sod or
	m_probEn == pda::Euler1d::Lax or
	m_probEn == pda::Euler1d::ShuOsher)
    {
      using ghost_filler_t  = pda::impl::Ghost1dNeumannFiller<
	U_t, MeshType, ghost_container_type>;
      ghost_filler_t ghF(stencilSizeNeeded, U, m_meshObj,
			 m_ghostLeft, m_ghostRight);

      const auto & rowsBd = m_meshObj.get().graphRowsOfCellsNearBd();
      for (decltype(rowsBd.size()) it=0; it<rowsBd.size(); ++it){
	ghF(rowsBd[it], it, m_numDofPerCell);
      }
    }
  }

  template<class U_t, class V_t>
  void velocityOnlyInnerCellsImpl(const U_t & U,
				  const scalar_type /*currentTime*/,
				  V_t & V,
				  flux_type & fluxL,
				  flux_type & fluxR,
				  edge_rec_type & uMinusHalfNeg,
				  edge_rec_type & uMinusHalfPos,
				  edge_rec_type & uPlusHalfNeg,
				  edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::ee::impl::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorForDiscreteFunction<
	    dimensionality, MeshType, U_t, edge_rec_type>,
	  3, scalar_type, flux_type>,
      V_t, scalar_type
      >;

    functor_type F(V, m_meshObj.get().dxInv(),
		   /* end args for velo */
		   m_inviscidFluxSchemeEn, m_gamma, fluxL, fluxR,
		   /* end args for flux */
		   toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj.get(),
		   uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		   /* end args for reconstructor */
		   );

    const auto & graphRows = m_meshObj.get().graphRowsOfCellsAwayFromBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      F(graphRows[it], m_numDofPerCell);
    }
  }

  template<class U_t, class V_t>
  void velocityAndJacobianImpl(const U_t & U,
			       const scalar_type currentTime,
			       V_t & V,
			       jacobian_type & J,
			       flux_type & fluxL,
			       flux_type & fluxR,
			       edge_rec_type & uMinusHalfNeg,
			       edge_rec_type & uMinusHalfPos,
			       edge_rec_type & uPlusHalfNeg,
			       edge_rec_type & uPlusHalfPos) const
  {
    // flux jacobians
    flux_jac_type fluxJacLNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacLPos(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacRNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacRPos(m_numDofPerCell, m_numDofPerCell);

    int nonZerosCountBeforeComputing = J.nonZeros();
    (void) nonZerosCountBeforeComputing;

    // this is done in parallel if openmp is on
    velocityAndJacInnerCellsImpl(U, currentTime, V, J,
				 fluxL, fluxR,
				 fluxJacLNeg, fluxJacLPos,
				 fluxJacRNeg, fluxJacRPos,
				 uMinusHalfNeg, uMinusHalfPos,
				 uPlusHalfNeg,  uPlusHalfPos);

    // near boundary I have be careful because
    // the jacobian can only be first order for now
    // also, in 1d there is not enough work to justify threading
    if (!m_meshObj.get().isFullyPeriodic())
      {
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp single
#endif
	{
	  if (m_inviscidFluxRecEn == InviscidFluxReconstruction::FirstOrder){
	    velocityAndJacNearBDCellsImplFirstOrder(U, currentTime, V, J,
						    fluxL, fluxR,
						    fluxJacLNeg, fluxJacLPos,
						    fluxJacRNeg, fluxJacRPos,
						    uMinusHalfNeg, uMinusHalfPos,
						    uPlusHalfNeg,  uPlusHalfPos);
	  }

	  else{
	    velocityAndJacNearBDCellsImplDifferentScheme(U, currentTime, V, J,
							 fluxL, fluxR,
							 fluxJacLNeg, fluxJacLPos,
							 fluxJacRNeg, fluxJacRPos,
							 uMinusHalfNeg, uMinusHalfPos,
							 uPlusHalfNeg,  uPlusHalfPos);
	  }
	}
    }

    assert(J.nonZeros() == nonZerosCountBeforeComputing);
  }


  template<class U_t, class V_t>
  void velocityAndJacInnerCellsImpl(const U_t & U,
				    const scalar_type /*currentTime*/,
				    V_t & V,
				    jacobian_type & J,
				    flux_type & fluxL,
				    flux_type & fluxR,
				    flux_jac_type & fluxJacLNeg,
				    flux_jac_type & fluxJacLPos,
				    flux_jac_type & fluxJacRNeg,
				    flux_jac_type & fluxJacRPos,
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

    // allocate gradients of reconstructed states
    // the size depends on the scheme selected
    const int stencilSize = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    reconstruction_gradient_t gradLNeg(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradLPos(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRNeg(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRPos(m_numDofPerCell, stencilSize-1);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impl::ComputeDirectionalFluxBalanceJacobianOnInteriorCell<
	  pda::ee::impl::ComputeDirectionalFluxValuesAndJacobians<
	    pda::impl::ReconstructorForDiscreteFunction<
	      dimensionality, MeshType, U_t, edge_rec_type, reconstruction_gradient_t>,
	    3, scalar_type, flux_type, flux_jac_type>,
	  dimensionality, MeshType, jacobian_type>,
      V_t, scalar_type
      >;

    functor_type F(V, m_meshObj.get().dxInv(),
		   /* end args for velo */
		   J, xAxis, m_meshObj.get(),
		   /* end args for jac */
		   m_inviscidFluxSchemeEn, m_gamma, fluxL, fluxR,
		   fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
		   /* end args for flux */
		   toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj.get(),
		   uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		   gradLNeg, gradLPos, gradRNeg, gradRPos
		   /* end args for reconstructor */
		   );

    const auto & graphRows = m_meshObj.get().graphRowsOfCellsAwayFromBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      F(graphRows[it], m_numDofPerCell);
    }
  }

  template<class U_t, class V_t>
  void velocityOnlyNearBdCellsImpl(const U_t & U,
				   const scalar_type /*currentTime*/,
				   V_t & V,
				   flux_type & fluxL,
				   flux_type & fluxR,
				   edge_rec_type & uMinusHalfNeg,
				   edge_rec_type & uMinusHalfPos,
				   edge_rec_type & uPlusHalfNeg,
				   edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;

    const auto stencilSize = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    stencil_container_type stencilVals(m_numDofPerCell*stencilSize);

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, stencil_container_type,
      U_t, MeshType, ghost_container_type>;
    stencil_filler_t FillStencilF(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				  U, m_meshObj.get(), m_ghostLeft, m_ghostRight,
				  stencilVals);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::ee::impl::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  3, scalar_type, flux_type>,
      V_t, scalar_type
      >;

    functor_type func(V, m_meshObj.get().dxInv(),
		      /* end args for velo */
		      m_inviscidFluxSchemeEn, m_gamma, fluxL, fluxR,
		      /* end args for flux */
		      toReconstructionScheme(m_inviscidFluxRecEn), stencilVals,
		      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		      /* end args for reconstructor */
		      );

    const auto & graphRows = m_meshObj.get().graphRowsOfCellsNearBd();
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      FillStencilF(smPt, it, m_numDofPerCell);
      func(smPt, m_numDofPerCell);
    }
  }

  template<class U_t, class V_t>
  void velocityAndJacNearBDCellsImplFirstOrder(const U_t & U,
					       const scalar_type /*currentTime*/,
					       V_t & V,
					       jacobian_type & J,
					       flux_type & fluxL,
					       flux_type & fluxR,
					       flux_jac_type & fluxJacLNeg,
					       flux_jac_type & fluxJacLPos,
					       flux_jac_type & fluxJacRNeg,
					       flux_jac_type & fluxJacRPos,
					       edge_rec_type & uMinusHalfNeg,
					       edge_rec_type & uMinusHalfPos,
					       edge_rec_type & uPlusHalfNeg,
					       edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    assert(m_inviscidFluxRecEn == InviscidFluxReconstruction::FirstOrder);

    // if here, then the scheme for velocity matches
    // the one for Jacobian so we can use same functors

    const auto stencilSize = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    stencil_container_type stencilVals(m_numDofPerCell*stencilSize);

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, stencil_container_type,
      U_t, MeshType, ghost_container_type>;
    stencil_filler_t FillStencilF(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				  U, m_meshObj.get(), m_ghostLeft, m_ghostRight,
				  stencilVals);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impl::ComputeDirectionalFluxBalanceFirstOrderJacobianOnBoundaryCell<
	  pda::ee::impl::ComputeDirectionalFluxValuesAndJacobians<
	    pda::impl::ReconstructorFromStencil<
	      edge_rec_type, stencil_container_type>,
	    3, scalar_type, flux_type, flux_jac_type>,
	  dimensionality, MeshType, jacobian_type>,
      V_t, scalar_type
      >;

    functor_type func(V, m_meshObj.get().dxInv(),
		      /* end args for velo */
		      J, xAxis, m_meshObj.get(),
		      /* end args for jac */
		      m_inviscidFluxSchemeEn, m_gamma, fluxL, fluxR,
		      fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
		      /* end args for flux */
		      toReconstructionScheme(m_inviscidFluxRecEn), stencilVals,
		      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		      /* end args for reconstructor */
		      );

    std::array<scalar_type, 3> bdCellJacFactors;
    // both Sod and Lax have Neumann BC type
    bdCellJacFactors.fill(static_cast<scalar_type>(1));

    const auto & graphRows = m_meshObj.get().graphRowsOfCellsNearBd();
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      FillStencilF(smPt, it, m_numDofPerCell);
      func(smPt, m_numDofPerCell, bdCellJacFactors);
    }
  }

  template<class U_t, class V_t>
  void velocityAndJacNearBDCellsImplDifferentScheme(const U_t & U,
						    const scalar_type /*currentTime*/,
						    V_t & V,
						    jacobian_type & J,
						    flux_type & fluxL,
						    flux_type & fluxR,
						    flux_jac_type & fluxJacLNeg,
						    flux_jac_type & fluxJacLPos,
						    flux_jac_type & fluxJacRNeg,
						    flux_jac_type & fluxJacRPos,
						    edge_rec_type & uMinusHalfNeg,
						    edge_rec_type & uMinusHalfPos,
						    edge_rec_type & uPlusHalfNeg,
						    edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;

    // if here, then the velocity must be computed with Weno,
    /// while the jacobian must be computed with first order

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, stencil_container_type,
      U_t, MeshType, ghost_container_type>;

    // *****************************
    // *** functors for velocity ***
    // *****************************

    const auto stencilSizeV = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    stencil_container_type stencilValsV(m_numDofPerCell*stencilSizeV);

    stencil_filler_t FillStencilVelo(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				     U, m_meshObj.get(), m_ghostLeft, m_ghostRight,
				     stencilValsV);

    using velo_functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::ee::impl::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  3, scalar_type, flux_type>,
      V_t, scalar_type
      >;

    velo_functor_type funcVelo(V, m_meshObj.get().dxInv(),
			       /* end args for velo */
			       m_inviscidFluxSchemeEn, m_gamma, fluxL, fluxR,
			       /* end args for flux */
			       toReconstructionScheme(m_inviscidFluxRecEn), stencilValsV,
			       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			       /* end args for reconstructor */
			       );

    // *****************************
    // *** functors for jacobian ***
    // *****************************
    const auto firstOrderRec = pda::InviscidFluxReconstruction::FirstOrder;
    const auto stencilSizeForJ = reconstructionTypeToStencilSize(firstOrderRec);
    stencil_container_type stencilValsForJ(m_numDofPerCell*stencilSizeForJ);
    stencil_filler_t FillStencilJac(stencilSizeForJ,
				    U, m_meshObj.get(), m_ghostLeft, m_ghostRight,
				    stencilValsForJ);

    using jac_functor_type =
      pda::impl::ComputeDirectionalFluxBalanceFirstOrderJacobianOnBoundaryCell<
	pda::ee::impl::ComputeDirectionalFluxJacobians<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  3, scalar_type, flux_jac_type>,
      dimensionality, MeshType, jacobian_type
      >;

    jac_functor_type funcJac(J, xAxis, m_meshObj.get(),
			     /* end args for jac */
			     m_inviscidFluxSchemeEn, m_gamma,
			     fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
			     /* end args for flux */
			     toReconstructionScheme(firstOrderRec), stencilValsForJ,
			     uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			     /* end args for reconstructor */
			     );

    // ************
    // loop
    // ************
    std::array<scalar_type, 3> bdCellJacFactors;
    bdCellJacFactors.fill(static_cast<scalar_type>(1));

    const auto & graphRows = m_meshObj.get().graphRowsOfCellsNearBd();
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      FillStencilVelo(smPt, it, m_numDofPerCell);
      funcVelo(smPt, m_numDofPerCell);
      FillStencilJac(smPt, it, m_numDofPerCell);
      funcJac(smPt, m_numDofPerCell, bdCellJacFactors);
    }
  }

  void allocateStencilValuesContainer()
  {
    // the stencil size needed is determined by the desired reconstruction
    // kind NOT from the mesh. THis is important because for example
    // the mesh can have a wider connectivity that what is needed.
    const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    ::pressiodemoapps::resize(m_stencilVals, m_numDofPerCell*stencilSizeNeeded);
  }

  void allocateGhostValues()
  {
    /*
      In 1d, the ghost values at left and right are stored in arrays.
      The size of this array depends on the stencil.
      For now, we has as many rows as BD cells so that for each
      BD cell, the corresponding ghosts can be easily found.
      This could be changed since in 1d we literally just need two arrays only.

      stencil = 3:
      ghost array has 3 entries for the 3 dofs in the ghost cell.
      For example at left boundary we have:
      ----------------||
      |	 0,  1,     2 ||
      | rho, rho*u, E ||
      |		      ||
      ----------------||

      stencil = 5:
      ghost array has 6 entries such that at left boundary we have:
      note the order of the indexing
      ----------------|---------------||
      |	 3,  4,     5 |	 0,  1,     2 ||
      | rho, rho*u, E | rho, rho*u, E ||
      |		      |		      ||
      ----------------|---------------||

      stencil = 7:
      ghost array has 9 entries such that at left boundary we have:
      note the order of the indexing
      ----------------|---------------|---------------||
      |	 6,  7,     8 |	 3,  4,     5 |	 0,  1,     2 ||
      | rho, rho*u, E | rho, rho*u, E | rho, rho*u, E ||
      |		      |		      |		      ||
      ----------------|---------------|---------------||
     */

    const auto stencilSize    = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    const auto numGhostValues = m_numDofPerCell*((stencilSize-1)/2);
    const index_t s1 = m_meshObj.get().numCellsBd();
    ::pressiodemoapps::resize(m_ghostLeft,  s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight, s1, numGhostValues);
  }

protected:
  // common to all problems
  int m_numDofPerCell = {};
  ::pressiodemoapps::Euler1d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_inviscidFluxRecEn;
  ::pressiodemoapps::InviscidFluxScheme m_inviscidFluxSchemeEn;
  scalar_type m_gamma = {};
  std::reference_wrapper<const MeshType> m_meshObj;

  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  mutable stencil_container_type m_stencilVals = {};
  mutable ghost_container_type   m_ghostLeft = {};
  mutable ghost_container_type   m_ghostRight = {};
};

template<class MeshType> constexpr int EigenApp<MeshType>::dimensionality;

} // end namespace impleuler1d
} // end namespace pressiodemoapps
#endif
