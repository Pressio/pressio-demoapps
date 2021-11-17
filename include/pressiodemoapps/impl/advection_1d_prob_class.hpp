
#ifndef PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_
#define PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_

#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "functor_cell_jacobian_first_order.hpp"
#include "functor_cell_jacobian_weno.hpp"

namespace pressiodemoapps{ namespace impladv{

template<typename sc_t>
void linAdvRusanovFlux(sc_t & F, const sc_t & qL, const sc_t & qR){
  // for linear advection it boils down to
  F = qL;
}

template<typename sc_t>
void linAdvRusanovFluxJacobian(sc_t & JL, sc_t & JR,
			       const sc_t & qL, const sc_t & qR)
{
  JL = static_cast<sc_t>(1);
  JR = static_cast<sc_t>(0);
}

template<class MeshType>
class EigenAdvection1dApp
{

public:
  using index_t	      = typename MeshType::index_t;
  using scalar_type   = typename MeshType::scalar_t;
  using state_type    = Eigen::Matrix<scalar_type,-1,1>;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, index_t>;

  static constexpr int	   dimensionality{1};
  static constexpr index_t numDofPerCell{1};

private:
  using stencil_values_t = Eigen::Matrix<scalar_type,-1,1>;
  using reconstruction_gradient_t = Eigen::Matrix<scalar_type,-1,1>;

public:
  EigenAdvection1dApp(const MeshType & meshObj,
		      ::pressiodemoapps::Advection1d probEnum,
		      ::pressiodemoapps::InviscidFluxReconstruction recEnum)
    : m_meshObj(meshObj),
      m_numDofStencilMesh(m_meshObj.stencilMeshSize()),
      m_numDofSampleMesh(m_meshObj.sampleMeshSize()),
      m_probEn(probEnum), m_recEn(recEnum)
  {
    const auto stencilSize = reconstructionTypeToStencilSize(recEnum);
    ::pressiodemoapps::resize(m_stencilVals, stencilSize);
  }

  state_type initialCondition() const{
    state_type res(m_numDofStencilMesh);
    const auto & x = m_meshObj.viewX();
    for (int i=0; i<::pressiodemoapps::extent(x,0); ++i){
      res(i) = std::sin(M_PI*x(i));
    }
    return res;
  }

protected:
  // note that here we MUST use a template because when doing
  // bindings, this gets deduced to be a Ref
  template<class U_t, class V_t>
  void velocityAndOptionalJacobian(const U_t & U,
				   const scalar_type currentTime,
				   V_t & V,
				   jacobian_type * J) const
  {

    constexpr int xAxis = 1;

    // reconstructions values
    scalar_type uMinusHalfNeg{0}, uMinusHalfPos{0};
    scalar_type uPlusHalfNeg {0}, uPlusHalfPos {0};
    // fluxes
    scalar_type fluxL{0}, fluxR{0};
    // flux jacobians
    scalar_type fluxJacLNeg, fluxJacLPos;
    scalar_type fluxJacRNeg, fluxJacRPos;

    // gradients of reconstructed states
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    reconstruction_gradient_t gradLNeg(stencilSize-1);
    reconstruction_gradient_t gradLPos(stencilSize-1);
    reconstruction_gradient_t gradRNeg(stencilSize-1);
    reconstruction_gradient_t gradRPos(stencilSize-1);

    if (m_recEn == InviscidFluxReconstruction::FirstOrder)
    {
      using reconstruct_functor_t = ::pressiodemoapps::impl::ReconstructorFromState<
	dimensionality, numDofPerCell, scalar_type, U_t, MeshType>;
      reconstruct_functor_t Reconstructor(xAxis, m_recEn, U, m_meshObj,
					  uMinusHalfNeg, uMinusHalfPos,
					  uPlusHalfNeg,  uPlusHalfPos);

      using jac_fnct_t = ::pressiodemoapps::impl::FirstOrderInnerCellJacobianFunctor<
	dimensionality, numDofPerCell, scalar_type, jacobian_type, scalar_type, MeshType>;
      jac_fnct_t CellJacobianFunctor(*J, m_meshObj, fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos);

      loopImpl(U, currentTime, V, J,
	       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos,
	       fluxL, fluxR, fluxJacLNeg, fluxJacLPos,
	       fluxJacRNeg, fluxJacRPos,
	       Reconstructor, CellJacobianFunctor);
    }

    else if (m_recEn == InviscidFluxReconstruction::Weno3 or
	     m_recEn == InviscidFluxReconstruction::Weno5)
    {

      if (J){
	using reconstruct_functor_t = ::pressiodemoapps::impl::ReconstructorFromState<
	  dimensionality, numDofPerCell, scalar_type, U_t, MeshType, reconstruction_gradient_t>;
	reconstruct_functor_t Reconstructor(xAxis, m_recEn, U, m_meshObj,
					    uMinusHalfNeg, uMinusHalfPos,
					    uPlusHalfNeg,  uPlusHalfPos,
					    gradLNeg, gradLPos,
					    gradRNeg, gradRPos);

	using jac_fnct_t = ::pressiodemoapps::impl::WenoInnerCellJacobianFunctor<
	  dimensionality, numDofPerCell, scalar_type, jacobian_type,
	  scalar_type, MeshType, reconstruction_gradient_t>;
	jac_fnct_t CellJacobianFunctor(m_recEn, *J, m_meshObj,
				       fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
				       gradLNeg, gradLPos, gradRNeg, gradRPos);

	loopImpl(U, currentTime, V, J,
		 uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos,
		 fluxL, fluxR, fluxJacLNeg, fluxJacLPos,
		 fluxJacRNeg, fluxJacRPos,
		 Reconstructor, CellJacobianFunctor);
      }
      else
      {
	using reconstruct_functor_t = ::pressiodemoapps::impl::ReconstructorFromState<
	  dimensionality, numDofPerCell, scalar_type, U_t, MeshType>;
	reconstruct_functor_t Reconstructor(xAxis, m_recEn, U, m_meshObj,
					    uMinusHalfNeg, uMinusHalfPos,
					    uPlusHalfNeg,  uPlusHalfPos);

	loopImpl(U, currentTime, V, J,
		 uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos,
		 fluxL, fluxR, fluxJacLNeg, fluxJacLPos,
		 fluxJacRNeg, fluxJacRPos,
		 Reconstructor,
		 // no-op lambda for cell jacobian since here the Jacobian is disabled
		 [](index_t /*unused*/){}
		 );
      }
    }
  }

  void initializeJacobian(jacobian_type & J)
  {
    J.resize(m_numDofSampleMesh, m_numDofStencilMesh);

    using Tr = Eigen::Triplet<scalar_type>;
    std::vector<Tr> trList;

    const scalar_type zero = 0;
    const auto & graph = m_meshObj.graph();
    for (int cell=0; cell<m_meshObj.sampleMeshSize(); ++cell)
      {
	const auto jacRowOfCurrentCell = cell*numDofPerCell;
	const auto L0 = graph(cell, 1);
	const auto R0 = graph(cell, 2);
	const auto ciL0 = L0*numDofPerCell;
	const auto ci  = graph(cell, 0)*numDofPerCell;
	const auto ciR0 = R0*numDofPerCell;

	trList.push_back( Tr(jacRowOfCurrentCell, ciL0, zero) );
	trList.push_back( Tr(jacRowOfCurrentCell, ci,   zero) );
	trList.push_back( Tr(jacRowOfCurrentCell, ciR0, zero) );

	if (m_recEn == InviscidFluxReconstruction::Weno3 or
	    m_recEn == InviscidFluxReconstruction::Weno5){
	  const auto L1 = graph(cell, 3);
	  const auto R1 = graph(cell, 4);
	  trList.push_back( Tr(jacRowOfCurrentCell, L1*numDofPerCell, zero) );
	  trList.push_back( Tr(jacRowOfCurrentCell, R1*numDofPerCell, zero) );
	}
	if (m_recEn == InviscidFluxReconstruction::Weno5){
	  const auto L2 = graph(cell, 5);
	  const auto R2 = graph(cell, 6);
	  trList.push_back( Tr(jacRowOfCurrentCell, L2*numDofPerCell, zero) );
	  trList.push_back( Tr(jacRowOfCurrentCell, R2*numDofPerCell, zero) );
	}
      }

    J.setFromTriplets(trList.begin(), trList.end());

    // compress to make it a real Crs matrix
    if (!J.isCompressed()){
      J.makeCompressed();
    }
  }

private:
  template<class U_t, class V_t, class ReconstructFunctor, class CellJacobianFunctor>
  void loopImpl(const U_t & U,
		const scalar_type currentTime,
		V_t & V,
		jacobian_type * J,
		scalar_type & uMinusHalfNeg,
		scalar_type & uMinusHalfPos,
		scalar_type & uPlusHalfNeg,
		scalar_type & uPlusHalfPos,
		scalar_type & fluxL,
		scalar_type & fluxR,
		scalar_type & fluxJacLNeg,
		scalar_type & fluxJacLPos,
		scalar_type & fluxJacRNeg,
		scalar_type & fluxJacRPos,
		ReconstructFunctor & reconstructor,
		CellJacobianFunctor && cellJacobianFunctor) const
  {
    int nonZerosCountBeforeComputing = 0;
    if (J){
      nonZerosCountBeforeComputing = J->nonZeros();
      ::pressiodemoapps::set_zero(*J);
    }

    // loop over cells
    const auto dxInv = m_meshObj.dxInv();
    const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
    {
      const auto vIndex = smPt*numDofPerCell;

      reconstructor(smPt);
      linAdvRusanovFlux(fluxL, uMinusHalfNeg, uMinusHalfPos);
      linAdvRusanovFlux(fluxR, uPlusHalfNeg,  uPlusHalfPos);

      if (J){
	linAdvRusanovFluxJacobian(fluxJacLNeg, fluxJacLPos, uMinusHalfNeg, uMinusHalfPos);
	linAdvRusanovFluxJacobian(fluxJacRNeg, fluxJacRPos, uPlusHalfNeg,  uPlusHalfPos);
      }

      V(vIndex) = dxInv*(fluxL - fluxR);

      if (J){
	cellJacobianFunctor(smPt);
      }
    }

#ifdef NDEBUG
    if (J){
      // ensure that nonzeros count does not change
      // this can happen for instance if during the Jacobian
      // evlauation, new non-zero elements get inserted
      // std::cout << "NONZEROS BEFORE COMP = " << nonZerosCountBeforeComputing << "\n";
      // std::cout << "NONZEROS AFTER  COMP = " << J->nonZeros() << "\n";
      assert(J->nonZeros() == nonZerosCountBeforeComputing);
    }
#endif
  }

protected:
  const MeshType & m_meshObj;
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};
  ::pressiodemoapps::Advection1d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  mutable stencil_values_t m_stencilVals;
};

}}//end namespace
#endif
