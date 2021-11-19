
#ifndef PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_
#define PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_

#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "advection_velocity_mixin.hpp"
#include "mixin_cell_jacobian.hpp"
#include "advection_flux_mixin.hpp"

namespace pressiodemoapps{ namespace impladv{

template<class MeshType>
class EigenAdvection1dApp
{

public:
  using index_t	      = typename MeshType::index_t;
  using scalar_type   = typename MeshType::scalar_t;
  using state_type    = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, index_t>;

  static constexpr int	   dimensionality{1};
  static constexpr index_t numDofPerCell{1};

private:
  using stencil_values_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using reconstruction_gradient_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;

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
    J.resize(m_numDofSampleMesh, m_numDofStencilMesh);

    using Tr = Eigen::Triplet<scalar_type>;
    std::vector<Tr> trList;

    const scalar_type zero = 0;
    const auto & graph = m_meshObj.graph();
    for (int cell=0; cell<m_meshObj.sampleMeshSize(); ++cell)
      {
	const auto jacRowOfCurrentCell = cell*numDofPerCell;
	const auto ci  = graph(cell, 0)*numDofPerCell;

	trList.push_back( Tr(jacRowOfCurrentCell, ci,   zero) );

	const auto numNeighbors =
	  (m_recEn == InviscidFluxReconstruction::FirstOrder) ? 2
	    : (m_recEn == InviscidFluxReconstruction::Weno3) ? 4 : 6;

	for (int i=1; i<=numNeighbors; ++i){
	  trList.push_back( Tr(jacRowOfCurrentCell, graph(cell, i), zero) );
	}

	// const auto L0 = graph(cell, 1);
	// const auto R0 = graph(cell, 2);
	// const auto ciL0 = L0*numDofPerCell;
	// const auto ciR0 = R0*numDofPerCell;
	// trList.push_back( Tr(jacRowOfCurrentCell, ciL0, zero) );
	// trList.push_back( Tr(jacRowOfCurrentCell, ciR0, zero) );
	// if (m_recEn == InviscidFluxReconstruction::Weno3 or
	//     m_recEn == InviscidFluxReconstruction::Weno5){
	//   const auto L1 = graph(cell, 3);
	//   const auto R1 = graph(cell, 4);
	//   trList.push_back( Tr(jacRowOfCurrentCell, L1*numDofPerCell, zero) );
	//   trList.push_back( Tr(jacRowOfCurrentCell, R1*numDofPerCell, zero) );
	// }
	// if (m_recEn == InviscidFluxReconstruction::Weno5){
	//   const auto L2 = graph(cell, 5);
	//   const auto R2 = graph(cell, 6);
	//   trList.push_back( Tr(jacRowOfCurrentCell, L2*numDofPerCell, zero) );
	//   trList.push_back( Tr(jacRowOfCurrentCell, R2*numDofPerCell, zero) );
	// }
      }

    J.setFromTriplets(trList.begin(), trList.end());

    // compress to make it a real Crs matrix
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
    V.setZero();
    int nonZerosCountBeforeComputing = 0;

    if (J){
      nonZerosCountBeforeComputing = J->nonZeros();
      ::pressiodemoapps::set_zero(*J);

      velocityAndJacImpl(U, currentTime, V, *J);

      // std::cout << "NONZEROS BEFORE COMP = " << nonZerosCountBeforeComputing << "\n";
      // std::cout << "NONZEROS AFTER  COMP = " << J->nonZeros() << "\n";
      assert(J->nonZeros() == nonZerosCountBeforeComputing);
    }
    else{
      velocityOnlyImpl(U, currentTime, V);
    }
  }

private:
  template<class U_t, class V_t>
  void velocityAndJacImpl(const U_t & U,
			  const scalar_type currentTime,
			  V_t & V,
			  jacobian_type & J) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;

    // reconstructions values
    scalar_type uMinusHalfNeg{0}, uMinusHalfPos{0};
    scalar_type uPlusHalfNeg {0}, uPlusHalfPos {0};
    // fluxes
    scalar_type fluxL{0}, fluxR{0};
    // flux jacobians
    scalar_type fluxJacLNeg, fluxJacLPos;
    scalar_type fluxJacRNeg, fluxJacRPos;

    // allocate gradients of reconstructed states
    // the size depends on the scheme selected
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    reconstruction_gradient_t gradLNeg(stencilSize-1);
    reconstruction_gradient_t gradLPos(stencilSize-1);
    reconstruction_gradient_t gradRNeg(stencilSize-1);
    reconstruction_gradient_t gradRPos(stencilSize-1);

    using functor_type =
      pda::impladv::CellVelocity<
	pda::impl::InnerCellJacobian<
	  pda::impladv::FluxValuesAndJacobians<
	    pda::impl::ReconstructorForDiscreteFunction<
	      dimensionality, numDofPerCell, MeshType, U_t, scalar_type, reconstruction_gradient_t>,
	    numDofPerCell, scalar_type, scalar_type, scalar_type>,
	  dimensionality, numDofPerCell, MeshType, jacobian_type>,
      V_t, scalar_type
      >;

    functor_type F(V, m_meshObj.dxInv(),
		   /* end args for velo */
		   J, xAxis, m_meshObj,
		   /* end args for jac */
		   m_fluxEn, fluxL, fluxR,
		   fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
		   /* end args for flux */
		   xAxis, toReconstructionScheme(m_recEn), U, m_meshObj,
		   uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		   gradLNeg, gradLPos, gradRNeg, gradRPos
		   /* end args for reconstructor */
		   );

    loopImpl(F);
  }

  template<class U_t, class V_t>
  void velocityOnlyImpl(const U_t & U,
			const scalar_type currentTime,
			V_t & V) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;

    // reconstructions values
    scalar_type uMinusHalfNeg{0}, uMinusHalfPos{0};
    scalar_type uPlusHalfNeg {0}, uPlusHalfPos {0};
    // fluxes
    scalar_type fluxL{0}, fluxR{0};
    // flux jacobians
    scalar_type fluxJacLNeg, fluxJacLPos;
    scalar_type fluxJacRNeg, fluxJacRPos;

    using functor_type =
      pda::impladv::CellVelocity<
	pda::impladv::FluxValues<
	  pda::impl::ReconstructorForDiscreteFunction<
	    dimensionality, numDofPerCell, MeshType, U_t, scalar_type>,
	  numDofPerCell, scalar_type, scalar_type>,
      V_t, scalar_type
      >;

    functor_type F(V, m_meshObj.dxInv(),
		   /* end args for velo */
		   m_fluxEn, fluxL, fluxR,
		   /* end args for flux */
		   xAxis, toReconstructionScheme(m_recEn), U, m_meshObj,
		   uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		   /* end args for reconstructor */
		   );
    loopImpl(F);
  }

  template<class F_t>
  void loopImpl(F_t & F) const
  {
    const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt){
      F(smPt);
    }
  }

protected:
  const MeshType & m_meshObj;
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};
  ::pressiodemoapps::Advection1d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  ::pressiodemoapps::InviscidFluxScheme m_fluxEn =
      ::pressiodemoapps::InviscidFluxScheme::Rusanov;

  mutable stencil_values_t m_stencilVals;
};

}}//end namespace
#endif
