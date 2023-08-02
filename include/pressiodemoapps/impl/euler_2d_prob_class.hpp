/*
//@HEADER
// ************************************************************************
//
// euler_2d_prob_class.hpp
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

#ifndef PRESSIODEMOAPPS_EULER2D_APP_HPP_
#define PRESSIODEMOAPPS_EULER2D_APP_HPP_

#include "noop.hpp"
#include "euler_rankine_hugoniot.hpp"
#include "euler_rusanov_flux_values_function.hpp"
#include "euler_rusanov_flux_jacobian_function.hpp"
#include "euler_2d_initial_condition.hpp"
#include "euler_2d_ghost_filler_neumann.hpp"
#include "euler_2d_ghost_filler_sedov2d_sym.hpp"
#include "euler_2d_ghost_filler_normal_shock.hpp"
#include "euler_2d_ghost_filler_double_mach_reflection.hpp"
#include "euler_2d_ghost_filler_cross_shock.hpp"
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

namespace pressiodemoapps{
namespace impleuler2d{

// tags are used inside he public create function: create_problem_...()
// to dispatch to the proper problem
// so add new ones if a new problem is added
struct TagCrossShock{};

template<
  class MeshType,
  class BCFunctorsHolderType = impl::NoOperation<void>
  >
class EigenApp
{

public:
  using index_t	      = typename MeshType::index_t;
  using scalar_type   = typename MeshType::scalar_t;
  using state_type    = Eigen::Matrix<scalar_type,-1,1>;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, index_t>;
  using mesh_connectivity_graph_type = typename MeshType::graph_t;

  static constexpr int dimensionality{2};
  static constexpr int numDofPerCell{4};

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
	   ::pressiodemoapps::Euler2d probEnum,
	   ::pressiodemoapps::InviscidFluxReconstruction recEnum,
	   ::pressiodemoapps::InviscidFluxScheme fluxEnum,
	   int icIdentifier)
    : m_meshObj(meshObj),
      m_probEn(probEnum),
      m_recEn(recEnum),
      m_fluxEn(fluxEnum),
      m_icIdentifier(icIdentifier)
  {
    m_numDofStencilMesh = m_meshObj.get().stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.get().sampleMeshSize() * numDofPerCell;
    allocateGhosts();
  }

#if !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  EigenApp(const MeshType & meshObj,
	   ::pressiodemoapps::Euler2d probEnum,
	   ::pressiodemoapps::InviscidFluxReconstruction recEnum,
	   ::pressiodemoapps::InviscidFluxScheme fluxEnum,
	   BCFunctorsHolderType && bcHolder,
	   int icIdentifier)
    : m_meshObj(meshObj),
      m_probEn(probEnum),
      m_recEn(recEnum),
      m_fluxEn(fluxEnum),
      m_icIdentifier(icIdentifier),
      m_bcFuncsHolder(std::move(bcHolder))
  {
    m_numDofStencilMesh = m_meshObj.get().stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.get().sampleMeshSize() * numDofPerCell;
    allocateGhosts();
  }
#endif

  EigenApp(TagCrossShock /*tag*/,
	   const MeshType & meshObj,
	   ::pressiodemoapps::InviscidFluxReconstruction recEnum,
	   ::pressiodemoapps::InviscidFluxScheme fluxEnum,
	   int icIdentifier,
	   scalar_type inletXVel,
	   scalar_type bottomYVel,
	   scalar_type density)
    : m_meshObj(meshObj),
      m_probEn(::pressiodemoapps::Euler2d::CrossShock),
      m_recEn(recEnum),
      m_fluxEn(fluxEnum),
      m_icIdentifier(icIdentifier),
      m_crossshock_params{density, inletXVel, bottomYVel}
  {
    m_numDofStencilMesh = m_meshObj.get().stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.get().sampleMeshSize() * numDofPerCell;
    allocateGhosts();
  }

  scalar_type gamma() const{
    return m_gamma;
  }

  state_type initialCondition() const{
    return initialConditionImpl();
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

    fillGhosts(U, currentTime);

    if (J){
      velocityAndJacobianImpl(U, currentTime, V, *J,
			      fluxL, fluxF, fluxR, fluxB,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);
    }

    else{
      if (!m_meshObj.get().isFullyPeriodic()){
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
	const auto jacRowOfCurrCellRho = smPt*numDofPerCell;

	// entries wrt current cell's dofs
	const auto jacColOfCurrCellRho = graph(smPt, 0)*numDofPerCell;
	for (int k=0; k<numDofPerCell; ++k){
	  for (int j=0; j<numDofPerCell; ++j){
	    trList.push_back( Tr(jacRowOfCurrCellRho+k, jacColOfCurrCellRho+j, zero) );
	  }
	}

	// wrt neighbors: this depends on the advection scheme
	const int numNeighbors =
	  (m_recEn == InviscidFluxReconstruction::FirstOrder) ? 4
	  : (m_recEn == InviscidFluxReconstruction::Weno3) ? 8
	  : (m_recEn == InviscidFluxReconstruction::Weno5) ? 12 : -1;
	assert(numNeighbors != -1);

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

  state_type initialConditionImpl() const
  {
    state_type initialState(m_numDofStencilMesh);

    switch(m_probEn)
      {
      case ::pressiodemoapps::Euler2d::PeriodicSmooth:{
	sin2dEulerIC(initialState, m_meshObj.get(), m_gamma);
	return initialState;
      }

      case ::pressiodemoapps::Euler2d::KelvinHelmholtz:{
	KelvinHelmholtzIC(initialState, m_meshObj.get(), m_gamma);
	return initialState;
      }

      case ::pressiodemoapps::Euler2d::SedovFull:{
	sedov2dIC(initialState, m_meshObj.get(), m_gamma);
	return initialState;
      }

      case ::pressiodemoapps::Euler2d::SedovSymmetry:{
	sedov2dsymmetryIC(initialState, m_meshObj.get(), m_gamma);
	return initialState;
      }

      case ::pressiodemoapps::Euler2d::Riemann:
#if !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
      case ::pressiodemoapps::Euler2d::RiemannCustomBCs:
#endif
      {
	if( m_icIdentifier == 1){
	  riemann2dIC1(initialState, m_meshObj.get(), m_gamma);
	  return initialState;
	}
	else if (m_icIdentifier == 2){
	  riemann2dIC2(initialState, m_meshObj.get(), m_gamma);
	  return initialState;
	}
	else{
	  throw std::runtime_error("Euler2d Riemann: invalid IC identifier");
	}
      }

      case ::pressiodemoapps::Euler2d::NormalShock:{
	normalShock2dIC(initialState, m_meshObj.get(), m_gamma);
	return initialState;
      }

      case ::pressiodemoapps::Euler2d::DoubleMachReflection:{
	doubleMachReflection2dIC(initialState, m_meshObj.get(), m_gamma);
	return initialState;
      }

      case ::pressiodemoapps::Euler2d::testingonlyneumann:{
	return initialState;
      }

      case ::pressiodemoapps::Euler2d::CrossShock:{
	crossShockIC(initialState, m_meshObj.get(), m_gamma,
		     m_crossshock_params[0], m_crossshock_params[1]);
	return initialState;
      }
      };

    return initialState;
  }


  template<class U_t>
  void fillGhosts(const U_t & U, const scalar_type currentTime) const
  {
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    if (m_probEn == ::pressiodemoapps::Euler2d::SedovFull or
	m_probEn == ::pressiodemoapps::Euler2d::Riemann or
	m_probEn == ::pressiodemoapps::Euler2d::testingonlyneumann)
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

#if !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
    else if (m_probEn == ::pressiodemoapps::Euler2d::RiemannCustomBCs)
    {
      const auto & x = m_meshObj.get().viewX();
      const auto & y = m_meshObj.get().viewY();
      const auto & graph = m_meshObj.get().graph();
      const auto & rowsBd = m_meshObj.get().graphRowsOfCellsNearBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
      for (decltype(rowsBd.size()) it=0; it<rowsBd.size(); ++it)
      {
	auto currentCellGraphRow = graph.row(rowsBd[it]);
	const int cellGID = rowsBd[it];
	const auto myX = x(cellGID);
	const auto myY = y(cellGID);

	/* IMPORTANT: keep the following as separate ifs wihtout ORs
	   because some cells might has ghosts on multiple sides so
	   we need these ifs not exclusive
	*/
	if (m_meshObj.get().hasBdLeft2d(cellGID)){
	  auto ghostVals = m_ghostLeft.row(it);
	  m_bcFuncsHolder(impl::GhostRelativeLocation::Left,
			  it, currentCellGraphRow, myX, myY, U, numDofPerCell,
			  m_meshObj.get().dx(), ghostVals);
	}

	if (m_meshObj.get().hasBdRight2d(cellGID)){
	  auto ghostVals = m_ghostRight.row(it);
	  m_bcFuncsHolder(impl::GhostRelativeLocation::Right,
			  it, currentCellGraphRow, myX, myY, U, numDofPerCell,
			  m_meshObj.get().dx(), ghostVals);
	}

	if (m_meshObj.get().hasBdBack2d(cellGID)){
	  auto ghostVals = m_ghostBack.row(it);
	  m_bcFuncsHolder(impl::GhostRelativeLocation::Back,
			  it, currentCellGraphRow, myX, myY, U, numDofPerCell,
			  m_meshObj.get().dy(), ghostVals);
	}

	if (m_meshObj.get().hasBdFront2d(cellGID)){
	  auto ghostVals = m_ghostFront.row(it);
	  m_bcFuncsHolder(impl::GhostRelativeLocation::Front,
			  it, currentCellGraphRow, myX, myY, U, numDofPerCell,
			  m_meshObj.get().dy(), ghostVals);
	}
      }
    }
#endif

    else if (m_probEn == ::pressiodemoapps::Euler2d::SedovSymmetry)
    {
      using ghost_filler_t  = Sedov2dSymmetryGhostFiller<
	U_t, MeshType, ghost_container_type>;
      ghost_filler_t ghF(stencilSize, U, m_meshObj.get(),
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

    else if (m_probEn == ::pressiodemoapps::Euler2d::NormalShock)
    {
      using ghost_filler_t = NormalShock2dGhostFiller<
	U_t, MeshType, ghost_container_type>;
      ghost_filler_t ghF(stencilSize, U,
			 currentTime, m_gamma, m_meshObj.get(),
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

    else if (m_probEn == ::pressiodemoapps::Euler2d::DoubleMachReflection)
    {
      using ghost_filler_t = DoubleMachReflection2dGhostFiller<
	U_t, MeshType, ghost_container_type>;
      ghost_filler_t ghF(stencilSize, U,
			 currentTime, m_gamma, m_meshObj.get(),
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

    else if (m_probEn == ::pressiodemoapps::Euler2d::CrossShock)
    {
      using ghost_filler_t  = CrossShock2dGhostFiller<
	U_t, MeshType, ghost_container_type>;
      ghost_filler_t ghF(stencilSize, U, m_gamma,
			 m_crossshock_params[0],
			 m_crossshock_params[1],
			 m_crossshock_params[2],
			 m_meshObj.get(),
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
    // only need to do near-BD cells if there are any
    if (!m_meshObj.get().isFullyPeriodic()){
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
    }

    velocityAndJacInnerCellsImpl(U, currentTime, V, J,
				 fluxL, fluxF, fluxR, fluxB,
				 fluxJacLNeg, fluxJacLPos,
				 fluxJacFNeg, fluxJacFPos,
				 fluxJacRNeg, fluxJacRPos,
				 fluxJacBNeg, fluxJacBPos,
				 uMinusHalfNeg, uMinusHalfPos,
				 uPlusHalfNeg,  uPlusHalfPos);

    // std::cout << J.nonZeros() << " "
    // 	      << nonZerosCountBeforeComputing << std::endl;
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

    const auto & graphRows = m_meshObj.get().graphRowsOfCellsAwayFromBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      Fx(smPt, numDofPerCell);
      Fy(smPt, numDofPerCell);
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
				  U, m_meshObj.get(), m_ghostLeft, m_ghostRight,
				  stencilVals, xAxis);

    stencil_filler_t FillStencilY(reconstructionTypeToStencilSize(m_recEn),
				  U, m_meshObj.get(), m_ghostBack, m_ghostFront,
				  stencilVals, yAxis);

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

    const auto & graphRows = m_meshObj.get().graphRowsOfCellsNearBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
      {
	const auto smPt = graphRows[it];

	FillStencilX(smPt, it, numDofPerCell);
	fillJacFactorsForCellBd(smPt, xAxis);
	funcx(smPt, numDofPerCell, m_bcCellJacFactors);

	FillStencilY(smPt, it, numDofPerCell);
	fillJacFactorsForCellBd(smPt, yAxis);
	funcy(smPt, numDofPerCell, m_bcCellJacFactors);
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
      dimensionality, stencil_container_type, U_t, MeshType, ghost_container_type>;

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

    // ************
    // loop
    // ************
    const auto & graphRows = m_meshObj.get().graphRowsOfCellsNearBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];

      FillStencilVeloX(smPt, it, numDofPerCell);
      funcVeloX(smPt, numDofPerCell);
      FillStencilJacX(smPt, it, numDofPerCell);
      fillJacFactorsForCellBd(smPt, xAxis);
      funcJacX(smPt, numDofPerCell, m_bcCellJacFactors);

      FillStencilVeloY(smPt, it, numDofPerCell);
      funcVeloY(smPt, numDofPerCell);
      FillStencilJacY(smPt, it, numDofPerCell);
      fillJacFactorsForCellBd(smPt, yAxis);
      funcJacY(smPt, numDofPerCell, m_bcCellJacFactors);
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

    const auto & graphRows = m_meshObj.get().graphRowsOfCellsAwayFromBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      Fx(smPt, numDofPerCell);
      Fy(smPt, numDofPerCell);
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
				  U, m_meshObj.get(), m_ghostLeft, m_ghostRight,
				  stencilVals, xAxis);

    stencil_filler_t FillStencilY(reconstructionTypeToStencilSize(m_recEn),
				  U, m_meshObj.get(), m_ghostBack, m_ghostFront,
				  stencilVals, yAxis);

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
    }
  }

  void fillJacFactorsForCellBd(index_t graphRow, int axis) const
  {
    assert(axis <= 2);

    if (m_probEn == ::pressiodemoapps::Euler2d::SedovFull or
	m_probEn == ::pressiodemoapps::Euler2d::Riemann or
	m_probEn == ::pressiodemoapps::Euler2d::testingonlyneumann)
    {
      (void)graphRow;
      // neumann
      m_bcCellJacFactors.fill(static_cast<scalar_type>(1));
      return;
    }

#if !defined PRESSIODEMOAPPS_ENABLE_BINDINGS
    else if (m_probEn == ::pressiodemoapps::Euler2d::RiemannCustomBCs)
    {
      const auto & x = m_meshObj.get().viewX();
      const auto & y = m_meshObj.get().viewY();
      const auto & graph = m_meshObj.get().graph();
      auto currentCellGraphRow = graph.row(graphRow);
      const int cellGID = currentCellGraphRow[0];
      const auto myX = x(cellGID);
      const auto myY = y(cellGID);

      if (axis==1){
	if (m_meshObj.get().hasBdLeft2d(graphRow)){
	  m_bcFuncsHolder(impl::GhostRelativeLocation::Left,
			  currentCellGraphRow, myX, myY, numDofPerCell,
			  /*axis, ,*/ m_bcCellJacFactors);
	}
	if (m_meshObj.get().hasBdRight2d(graphRow)){
	  m_bcFuncsHolder(impl::GhostRelativeLocation::Right,
			  currentCellGraphRow, myX, myY, numDofPerCell,
			  /*axis, ,*/ m_bcCellJacFactors);
	}
      }
      else{
	if (m_meshObj.get().hasBdBack2d(graphRow)){
	  m_bcFuncsHolder(impl::GhostRelativeLocation::Back,
			  currentCellGraphRow, myX, myY, numDofPerCell,
			  /*axis, ,*/ m_bcCellJacFactors);
	}

	if (m_meshObj.get().hasBdFront2d(graphRow)){
	  m_bcFuncsHolder(impl::GhostRelativeLocation::Front,
			  currentCellGraphRow, myX, myY, numDofPerCell,
			  /*axis, ,*/ m_bcCellJacFactors);
	}
      }

      return;
    }
#endif

    else if (m_probEn == ::pressiodemoapps::Euler2d::SedovSymmetry)
    {
      if (axis == 1 && m_meshObj.get().hasBdLeft2d(graphRow)){
	// relfective
	m_bcCellJacFactors = {1., -1., 1., 1.};
	return;
      }

      else if (axis == 2 && m_meshObj.get().hasBdBack2d(graphRow)){
	// relfective
	m_bcCellJacFactors = {1., 1., -1., 1.};
	return;
      }

      else{
	// neumann
	m_bcCellJacFactors = {1., 1., 1., 1.};
	return;
      }
    }

    else if (m_probEn == ::pressiodemoapps::Euler2d::NormalShock)
    {
      if (axis == 2){
	// relfective
	m_bcCellJacFactors = {1., 1., -1., 1.};
	return;
      }
      else{
	// neumann
	m_bcCellJacFactors = {1., 1., 1., 1.};
	return;
      }
    }

    else if (m_probEn == ::pressiodemoapps::Euler2d::DoubleMachReflection)
    {

      const scalar_type wedgePosition = static_cast<scalar_type>(1)/static_cast<scalar_type>(6);
      const auto & x = m_meshObj.get().viewX();
      const auto cellGID = m_meshObj.get().graph()(graphRow, 0);
      const auto myX = x(cellGID);

      if (axis == 1){
	// homog neumann
	m_bcCellJacFactors = {1., 1., 1., 1.};
	return;
      }

      if (axis == 2 && m_meshObj.get().hasBdBack2d(graphRow) && (myX < wedgePosition)){
	// homog neumann
	m_bcCellJacFactors = {1., 1., 1., 1.};
	return;
      }

      if (axis == 2 && m_meshObj.get().hasBdBack2d(graphRow) && (myX >= wedgePosition)){
	// reflective
	m_bcCellJacFactors = {1., 1., -1., 1.};
	return;
      }

      // dirichlet
      m_bcCellJacFactors = {0., 0., 0., 0.};
      return;
    }

    else if (m_probEn == ::pressiodemoapps::Euler2d::CrossShock)
    {
      if (axis == 1 && m_meshObj.get().hasBdLeft2d(graphRow)){
	// dirichlet
	m_bcCellJacFactors = {0., 0., 0., 0.};
	return;
      }

      if (axis == 1 && m_meshObj.get().hasBdRight2d(graphRow)){
	// homog neumann
	m_bcCellJacFactors = {1., 1., 1., 1.};
	return;
      }

      if (axis == 2 && m_meshObj.get().hasBdBack2d(graphRow))
      {
	// dirichlet for u,v, and homog neumann for rho and rho*E
	m_bcCellJacFactors = {1., 0., 0., 1.};
	return;
      }

      if (axis == 2 && m_meshObj.get().hasBdFront2d(graphRow))
      {
	// homog neumann for rho and rho*E
	m_bcCellJacFactors = {1., 1., 1., 1.};
	return;
      }

    }
  }

  void allocateGhosts()
  {
    /*
      stencil = 3:
      ghost array has entries for the dofs in the ghost cell.
      For example at left boundary we have:
      ------------------------||
      |	 0,  1,     2,     3  ||
      | rho, rho*u, rho*v, E  ||
      |			      ||
      ------------------------||

      stencil = 5:
      ghost array has entries such that at left boundary we have:
      note the order of the indexing
      -----------------------|-----------------------||
      |	 4,  5,     6,     7 |	 0,  1,     2,     3 ||
      | rho, rho*u, rho*v, E | rho, rho*u, rho*v,  E ||
      |		             |		             ||
      -----------------------|-----------------------||

      stencil = 7:
      ghost array has 9 entries such that at left boundary we have:
      note the order of the indexing
      ------------------------|----------------------|----------------------||
      |	 8,    9,    10,   11 |	 4,    5,     6,   7 |	0,    1,     2,   3 ||
      | rho, rho*u, rho*v, E  | rho, rho*u, rho*v, E | rho, rho*u, rho*v, E ||
      |		              |		             |		            ||
      ------------------------|----------------------|----------------------||
     */

    const auto stencilSize    = reconstructionTypeToStencilSize(m_recEn);
    const auto numGhostValues = numDofPerCell*((stencilSize-1)/2);

    const index_t s1 = m_meshObj.get().numCellsBd();
    ::pressiodemoapps::resize(m_ghostLeft,  s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostFront, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostBack,  s1, numGhostValues);
  }

protected:
  scalar_type m_gamma = static_cast<scalar_type>(1.4);

  std::reference_wrapper<const MeshType> m_meshObj;
  ::pressiodemoapps::Euler2d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  ::pressiodemoapps::InviscidFluxScheme m_fluxEn;

  // which initial condition to use by default
  int m_icIdentifier = 1;

  // note that dof refers to the degress of freedom,
  // which is NOT same as grid points.
  // SampleMesh_ identifies the velocity/residual locations
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  mutable ghost_container_type m_ghostLeft;
  mutable ghost_container_type m_ghostFront;
  mutable ghost_container_type m_ghostRight;
  mutable ghost_container_type m_ghostBack;

  std::array<scalar_type, 2> normalX_{1, 0};
  std::array<scalar_type, 2> normalY_{0, 1};

  // for cross-shock problem: density, inletXVel, bottomYVel
  std::array<scalar_type, 3> m_crossshock_params;

  mutable std::array<scalar_type, numDofPerCell> m_bcCellJacFactors;

  BCFunctorsHolderType m_bcFuncsHolder = {};
};

template<class T1, class T2> constexpr int EigenApp<T1,T2>::numDofPerCell;
template<class T1, class T2> constexpr int EigenApp<T1,T2>::dimensionality;

}}//end namespace
#endif
