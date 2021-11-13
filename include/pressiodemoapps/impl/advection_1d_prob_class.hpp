
#ifndef PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_
#define PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_

#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "cell_jacobian_first_order.hpp"

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
  using stencil_values_t = state_type;

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
  void initializeJacobian(jacobian_type & J)
  {
    initializeJacobianFirstOrder(J);
    // compress to make it a real Crs matrix
    if (!J.isCompressed()){
      J.makeCompressed();
    }
  }

  void initializeJacobianFirstOrder(jacobian_type & J)
  {
    J.resize(m_numDofSampleMesh, m_numDofStencilMesh);

    using Tr = Eigen::Triplet<scalar_type>;
    std::vector<Tr> trList;

    const scalar_type val0 = 0;
    const auto & graph = m_meshObj.graph();
    for (int cell=0; cell<m_meshObj.sampleMeshSize(); ++cell)
      {
	const auto jacRowOfCurrentCell = cell*numDofPerCell;
	const auto L0 = graph(cell, 1);
	const auto R0 = graph(cell, 2);
	const auto ciL0 = L0*numDofPerCell;
	const auto ci  = graph(cell, 0)*numDofPerCell;
	const auto ciR0 = R0*numDofPerCell;

	trList.push_back( Tr(jacRowOfCurrentCell, ciL0, val0) );
	trList.push_back( Tr(jacRowOfCurrentCell, ci,   val0) );
	trList.push_back( Tr(jacRowOfCurrentCell, ciR0, val0) );
      }

    J.setFromTriplets(trList.begin(), trList.end());
  }

  // note that here we MUST use a template because when doing
  // bindings, this gets deduced to be a Ref
  template<class U_t, class V_t>
  void velocityAndOptionalJacobian(const U_t & U,
				   const scalar_type currentTime,
				   V_t & V,
				   jacobian_type * J) const
  {
    int nonZerosCountBeforeComputing = 0;
    if (J){
      nonZerosCountBeforeComputing = J->nonZeros();
      ::pressiodemoapps::set_zero(*J);
    }

    scalar_type FL{0}, FR{0};
    scalar_type uMinusHalfNeg{0}, uMinusHalfPos{0};
    scalar_type uPlusHalfNeg {0}, uPlusHalfPos {0};
    scalar_type JLneg, JLpos, JRneg, JRpos;

    // reconstruct functor for face fluxes
    // here we need to use whatever order (m_recEn) user decides
    using reconstruct_functor_t = ::pressiodemoapps::impl::ReconstructorFromState<
      dimensionality, scalar_type, U_t, MeshType>;
    reconstruct_functor_t Reconstructor(m_recEn, U, m_meshObj,
					uMinusHalfNeg, uMinusHalfPos,
					uPlusHalfNeg,  uPlusHalfPos);

    // cell jacobian functor
    // currently, REGARDLESS of the reconstruction scheme,
    // we have only first-order Jacobian so we need a first-order reconstructor
    // for the Jacobian
    const auto firstOrder = ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder;
    scalar_type uMinusHalfNegForJ;
    scalar_type uMinusHalfPosForJ;
    scalar_type uPlusHalfNegForJ;
    scalar_type uPlusHalfPosForJ;
    reconstruct_functor_t ReconstructorForJ(firstOrder, U, m_meshObj,
					    uMinusHalfNegForJ, uMinusHalfPosForJ,
					    uPlusHalfNegForJ,  uPlusHalfPosForJ);

    // for this problem, due to periodic BC, any cell is an "inner" cell
    using jac_fnct_t = impl::FirstOrderInnerCellJacobianFunctor<
      dimensionality, numDofPerCell, scalar_type, jacobian_type, scalar_type, MeshType>;
    jac_fnct_t CellJacobianFunctor(*J, m_meshObj, JLneg, JLpos, JRneg, JRpos);

    // loop over cells
    const auto dxInv = m_meshObj.dxInv();
    const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
    {
      const auto vIndex = smPt*numDofPerCell;

      Reconstructor.template operator()<numDofPerCell>(smPt);

      if (J){
	// for now, REGARDLESS of the reconstruction scheme,
	// we currently only have only first-order Jacobian so we need
	// to run the reconstructor for the Jacobian
	// which will ensure that uMinusNegForJ, etc have the right values
	ReconstructorForJ.template operator()<numDofPerCell>(smPt);
      }

      linAdvRusanovFlux(FL, uMinusHalfNeg, uMinusHalfPos);
      linAdvRusanovFlux(FR, uPlusHalfNeg,  uPlusHalfPos);

      if (J){
	linAdvRusanovFluxJacobian(JLneg, JLpos, uMinusHalfNegForJ, uMinusHalfPosForJ);
	linAdvRusanovFluxJacobian(JRneg, JRpos, uPlusHalfNegForJ,  uPlusHalfPosForJ);
      }

      V(vIndex) = dxInv*(FL - FR);

      if (J){
	CellJacobianFunctor(smPt);
      }
    }

#ifdef NDEBUG
    if (J){
      // ensure that nonzeros count does not change
      // this can happen for instance if during the Jacobian
      // evlauation, new non-zero elements get inserted
      assert(J->nonZeros() != nonZerosCountBeforeComputing);
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
