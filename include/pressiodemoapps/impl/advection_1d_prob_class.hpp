/*
//@HEADER
// ************************************************************************
//
// advection_1d_prob_class.hpp
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

#ifndef PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_
#define PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_

#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "advection_1d_mixins.hpp"
#include "mixin_directional_flux_balance.hpp"
#include "mixin_directional_flux_balance_jacobian.hpp"

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#include <omp.h>
#endif

namespace pressiodemoapps{
namespace impladvection1d{

// tags are used inside he public create function: create_problem_...()
// to dispatch to the proper problem
// so add new ones if a new problem is added
struct TagLinearAdvection{};

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
  using flux_type	          = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using edge_rec_type	          = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using flux_jac_type             = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;
  using reconstruction_gradient_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  EigenApp() = delete;

  EigenApp(TagLinearAdvection /*tag*/,
	   const MeshType & meshObj,
	   ::pressiodemoapps::InviscidFluxReconstruction inviscidFluxRecEn,
	   ::pressiodemoapps::InviscidFluxScheme invFluxSchemeEn,
	   scalar_type velocity,
	   int icIdentifier)
    : m_numDofPerCell(1),
      m_probEn(::pressiodemoapps::Advection1d::PeriodicLinear),
      m_inviscidFluxRecEn(inviscidFluxRecEn),
      m_inviscidFluxSchemeEn(invFluxSchemeEn),
      m_meshObj(meshObj),
      m_icIdentifier(icIdentifier),
      m_linear_adv_vel(velocity)
  {
    m_numDofStencilMesh = m_meshObj.get().stencilMeshSize() * m_numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.get().sampleMeshSize() * m_numDofPerCell;
  }

  state_type initialCondition() const{
    state_type res(m_numDofStencilMesh);
    const auto & x = m_meshObj.get().viewX();

    for (int i=0; i<::pressiodemoapps::extent(x,0); ++i){
      if (m_icIdentifier == 1){
	res(i) = std::sin(M_PI*x(i));
      }

      else if (m_icIdentifier == 2){
	const auto x1= static_cast<scalar_type>(1.2);
	const auto A1= static_cast<scalar_type>(200.0);
	const auto delta1 = 16.0;
	const auto x2= static_cast<scalar_type>(2.5);
	const auto A2= static_cast<scalar_type>(100.0);
	const auto delta2 = 36.0;

	const auto dx1Sq = (x(i)-x1)*(x(i)-x1);
	const auto dx2Sq = (x(i)-x2)*(x(i)-x2);
        res(i) = 0.8*std::exp(-A1*dx1Sq/delta1) + std::exp(-A2*dx2Sq/delta2);
      }

      else if (m_icIdentifier == 3){
	const auto x1= static_cast<scalar_type>(2.);
	const auto x2= static_cast<scalar_type>(3.0);
	const auto delta = 0.5*0.5;
	const auto dx1Sq = (x(i)-x1)*(x(i)-x1);
	const auto dx2Sq = (x(i)-x2)*(x(i)-x2);
        res(i) = std::exp(-dx1Sq/delta) + 0.5*std::exp(-dx2Sq/delta);
      }

      else if (m_icIdentifier == 4){
	const auto x1= static_cast<scalar_type>(1);
	const auto x2= static_cast<scalar_type>(3);
	const auto sf= static_cast<scalar_type>(8);
	const auto dx1 = (x(i)-x1);
	const auto dx2 = (x(i)-x2);
        res(i) = std::tanh(sf*dx1) - std::tanh(sf*dx2);
      }

      else{
	throw std::runtime_error("advection1d: invalid ic");
      }
    }
    return res;
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

    const auto zero = static_cast<scalar_type>(0);
    const auto & graph = m_meshObj.get().graph();
    for (int cell=0; cell<m_meshObj.get().sampleMeshSize(); ++cell)
      {
	const auto jacRowOfCurrentCell = cell*m_numDofPerCell;
	const auto ci  = graph(cell, 0)*m_numDofPerCell;

	// entry wrt current cell
	trList.push_back( Tr(jacRowOfCurrentCell, ci, zero) );

	const int numNeighbors =
	  (m_inviscidFluxRecEn == InviscidFluxReconstruction::FirstOrder) ? 2
	   : (m_inviscidFluxRecEn == InviscidFluxReconstruction::Weno3) ? 4
	    : (m_inviscidFluxRecEn == InviscidFluxReconstruction::Weno5) ? 6 : -1;
	assert(numNeighbors != -1);
	for (int i=1; i<=numNeighbors; ++i){
	  trList.push_back( Tr(jacRowOfCurrentCell, graph(cell, i), zero) );
	}
      }

    J.setFromTriplets(trList.begin(), trList.end());

    // compress to make it a real Crs matrix
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
      int nonZerosCountBeforeComputing = J->nonZeros();
      velocityAndJacImpl(U, currentTime, V, *J);
      assert(J->nonZeros() == nonZerosCountBeforeComputing);
      (void) nonZerosCountBeforeComputing;
    }
    else{
      velocityOnlyImpl(U, currentTime, V);
    }

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
} // end pragma omp
#endif
  }

private:
  template<class U_t, class V_t>
  void velocityAndJacImpl(const U_t & U,
			  const scalar_type currentTime,
			  V_t & V,
			  jacobian_type & J) const
  {
    (void) currentTime;

    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;

    // for omp, these are private variables for each thread
    // edge reconstructions
    edge_rec_type uMinusHalfNeg(m_numDofPerCell);
    edge_rec_type uMinusHalfPos(m_numDofPerCell);
    edge_rec_type uPlusHalfNeg(m_numDofPerCell);
    edge_rec_type uPlusHalfPos(m_numDofPerCell);
    // fluxes
    flux_type fluxL(m_numDofPerCell);
    flux_type fluxR(m_numDofPerCell);

    // flux jacobians
    flux_jac_type fluxJacLNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacLPos(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacRNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacRPos(m_numDofPerCell, m_numDofPerCell);

    // allocate gradients of reconstructed states
    // the size depends on the scheme selected
    const auto stencilSize = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    reconstruction_gradient_t gradLNeg(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradLPos(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRNeg(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRPos(m_numDofPerCell, stencilSize-1);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impl::ComputeDirectionalFluxBalanceJacobianOnInteriorCell<
	  pda::impladvection1d::ComputeDirectionalFluxValuesAndJacobians<
	    pda::impl::ReconstructorForDiscreteFunction<
	      dimensionality, MeshType, U_t, edge_rec_type, reconstruction_gradient_t>,
	    scalar_type, flux_type, flux_jac_type>,
	  dimensionality, MeshType, jacobian_type>,
      V_t, scalar_type
      >;

    functor_type F(V, m_meshObj.get().dxInv(),
		   /* end args for velo */
		   J, xAxis, m_meshObj.get(),
		   /* end args for jac */
		   m_inviscidFluxSchemeEn, fluxL, fluxR,
		   fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos, m_linear_adv_vel,
		   /* end args for flux */
		   toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj.get(),
		   uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		   gradLNeg, gradLPos, gradRNeg, gradRPos
		   /* end args for reconstructor */
		   );

    const auto sampleMeshSize = m_meshObj.get().sampleMeshSize();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt){
      F(smPt, m_numDofPerCell);
    }
  }

  template<class U_t, class V_t>
  void velocityOnlyImpl(const U_t & U,
			const scalar_type /*currentTime*/,
			V_t & V) const
  {
    namespace pda = ::pressiodemoapps;

    // for omp, these are private variables for each thread
    // edge reconstructions
    edge_rec_type uMinusHalfNeg(m_numDofPerCell);
    edge_rec_type uMinusHalfPos(m_numDofPerCell);
    edge_rec_type uPlusHalfNeg(m_numDofPerCell);
    edge_rec_type uPlusHalfPos(m_numDofPerCell);
    // fluxes
    flux_type fluxL(m_numDofPerCell);
    flux_type fluxR(m_numDofPerCell);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impladvection1d::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorForDiscreteFunction<
	    dimensionality, MeshType, U_t, edge_rec_type>,
	  scalar_type, flux_type>,
      V_t, scalar_type
      >;

    functor_type F(V, m_meshObj.get().dxInv(),
		   /* end args for velo */
		   m_inviscidFluxSchemeEn, fluxL, fluxR, m_linear_adv_vel,
		   /* end args for flux */
		   toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj.get(),
		   uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		   /* end args for reconstructor */
		   );

    const auto sampleMeshSize = m_meshObj.get().sampleMeshSize();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt){
      F(smPt, m_numDofPerCell);
    }
  }

protected:
  // common to all problems
  int m_numDofPerCell = {};
  ::pressiodemoapps::Advection1d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_inviscidFluxRecEn;
  ::pressiodemoapps::InviscidFluxScheme m_inviscidFluxSchemeEn;
  std::reference_wrapper<const MeshType> m_meshObj;
  int m_icIdentifier = {};

  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  // parameters specific to problems
  // will need to handle this better later
  scalar_type m_linear_adv_vel = {};
};

template<class MeshType> constexpr int EigenApp<MeshType>::dimensionality;

}}//end namespace
#endif
