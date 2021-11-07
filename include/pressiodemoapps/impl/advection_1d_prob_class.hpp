
#ifndef PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_
#define PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_

// note that the code below is pretty ugly for now, but it works.
// this is inside impl namespace for a reason, and will need to be
// improved later on but we have a starting point.

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

template<
  class ScalarType,
  class MeshType,
  class StateType,
  class VelocityType
  >
class Advection1dAppRhsOnly
{
public:
  using index_t		 = typename MeshType::index_t;
  using scalar_type	 = ScalarType;
  using state_type	 = StateType;
  using velocity_type	 = VelocityType;
  using stencil_values_t = StateType;

  static constexpr int dimensionality{1};
  static constexpr int numDofPerCell{1};

public:
  Advection1dAppRhsOnly(const MeshType & meshObj,
		      ::pressiodemoapps::Advection1d probEnum,
		      ::pressiodemoapps::InviscidFluxReconstruction recEnum)
    : m_meshObj(meshObj), m_probEn(probEnum), m_recEn(recEnum)
  {
    m_numDofStencilMesh = m_meshObj.stencilMeshSize();
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize();

    const auto stencilSize = reconstructionTypeToStencilSize(recEnum);
    ::pressiodemoapps::resize(m_stencilVals, stencilSize);
  }

  index_t totalDofSampleMesh()  const{ return m_numDofSampleMesh;  }
  index_t totalDofStencilMesh() const{ return m_numDofStencilMesh; }

  state_type initialCondition() const{
    return initialConditionImpl();
  }

  velocity_type createVelocity() const {
    velocity_type V(m_numDofSampleMesh);
    return V;
  }

  void velocity(const state_type & stateO,
		const scalar_type timeValue,
		velocity_type & veloO) const{
    velocityImpl(stateO, timeValue, veloO);
  }

private:
  state_type initialConditionImpl() const{
    state_type res(m_numDofStencilMesh);
    const auto & x = m_meshObj.viewX();
    for (int i=0; i<::pressiodemoapps::extent(x,0); ++i){
      res(i) = std::sin(M_PI*x(i));
    }
    return res;
  }

  void velocityImpl(const state_type & U,
		    const scalar_type t,
		    velocity_type & V) const
  {
    scalar_type FL{0};
    scalar_type FR{0};
    scalar_type uMinusHalfNeg{0};
    scalar_type uMinusHalfPos{0};
    scalar_type uPlusHalfNeg {0};
    scalar_type uPlusHalfPos {0};

    using reconstruct_functor_t = ::pressiodemoapps::impl::ReconstructorFromState<
      dimensionality, scalar_type, state_type, MeshType>;

    reconstruct_functor_t Reconstructor(m_recEn, U, m_meshObj,
					uMinusHalfNeg, uMinusHalfPos,
					uPlusHalfNeg,  uPlusHalfPos);

    const auto dxInv = m_meshObj.dxInv();
    for (index_t smPt=0; smPt < m_meshObj.sampleMeshSize(); ++smPt)
    {
      Reconstructor.template operator()<numDofPerCell>(smPt);
      linAdvRusanovFlux(FL, uMinusHalfNeg, uMinusHalfPos);
      linAdvRusanovFlux(FR, uPlusHalfNeg,  uPlusHalfPos);
      const auto vIndex = smPt*numDofPerCell;
      V(vIndex) = dxInv*(FL - FR);
    }
  }

protected:
  const MeshType & m_meshObj;
  ::pressiodemoapps::Advection1d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};
  mutable stencil_values_t m_stencilVals;
};


#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class ScalarType, class MeshType>
using EigenAdvection1dAppRhsOnly = Advection1dAppRhsOnly<ScalarType, MeshType,
							 Eigen::Matrix<ScalarType,-1,1>,
							 Eigen::Matrix<ScalarType,-1,1>
							 >;

template<class ScalarType, class MeshType>
class EigenAdvection1dAppWithJacobian
  : public EigenAdvection1dAppRhsOnly<ScalarType, MeshType>
{
  using base_t = EigenAdvection1dAppRhsOnly<ScalarType, MeshType>;

public:
  using typename base_t::index_t;
  using typename base_t::scalar_type;
  using typename base_t::state_type;
  using typename base_t::velocity_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, int>;

  static constexpr int	   dimensionality{1};
  static constexpr index_t numDofPerCell{1};

private:
  using typename base_t::stencil_values_t;

public:
  template<class ...Args>
  EigenAdvection1dAppWithJacobian(Args && ... args)
    : base_t(std::forward<Args>(args)...),
      m_jacobian(m_numDofSampleMesh, m_numDofStencilMesh)
  {
    initializeJacobian();
  }

  // the Jacobian is by default fused with the velocity,
  // this method allows one to disable the jacobian
  // so only velocity is computed
  void disableJacobian() {
    m_onlyComputeVelocity = true;
  }

  jacobian_type createJacobian() const {
    return m_jacobian;
  }

  void velocity(const state_type & state,
		const scalar_type currentTime,
		velocity_type & V) const
  {
    if (m_onlyComputeVelocity){
      base_t::velocity(state, currentTime, V);
    }
    else{

      zeroOutJacobianEntries();
      velocityAndJacImpl(state, currentTime, V);

#ifdef NDEBUG
      // ensure that nonzeros count does not change
      // this can happen for instance if during the Jacobian
      // evlauation, new non-zero elements get inserted
      assert(m_jacNonZerosCount == m_jacobian.nonZeros());
#endif
    }
  }

  void jacobian(const state_type & state,
		const scalar_type timeValue,
		jacobian_type & J) const
  {
    if (!m_onlyComputeVelocity){
      // relies on jacobian been computed in velocity
      J = m_jacobian;
    }
  }

private:
  void zeroOutJacobianEntries() const{
    auto values = m_jacobian.valuePtr();
    for (int i=0; i<m_jacobian.nonZeros(); ++i){
      values[i] = 0.0;
    }
  }

  void initializeJacobian()
  {
    initializeJacobianFirstOrder();

    // compress to make it a real Crs matrix
    if (!m_jacobian.isCompressed()){
      m_jacobian.makeCompressed();
    }

    m_jacNonZerosCount = m_jacobian.nonZeros();

    // if Jacobian is disabled, free it
    if (m_onlyComputeVelocity){
      ::pressiodemoapps::resize(m_jacobian, 0, 0);
    }
  }

  void initializeJacobianFirstOrder()
  {
    using Tr = Eigen::Triplet<ScalarType>;
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
	trList.push_back( Tr(jacRowOfCurrentCell, ci, val0) );
	trList.push_back( Tr(jacRowOfCurrentCell, ciR0, val0) );
      }

    m_jacobian.setFromTriplets(trList.begin(), trList.end());
  }

  void velocityAndJacImpl(const state_type & U,
			  const scalar_type currentTime,
			  velocity_type & V) const
  {

    scalar_type FL{0};
    scalar_type FR{0};
    scalar_type uMinusHalfNeg{0};
    scalar_type uMinusHalfPos{0};
    scalar_type uPlusHalfNeg {0};
    scalar_type uPlusHalfPos {0};

    scalar_type JLneg;
    scalar_type JLpos;
    scalar_type JRneg;
    scalar_type JRpos;

    // reconstruct functor for face fluxes
    // here we need to use whatever order (m_recEn) user decides
    using reconstruct_functor_t = ::pressiodemoapps::impl::ReconstructorFromState<
      base_t::dimensionality, scalar_type, state_type, MeshType>;
    reconstruct_functor_t Reconstructor(m_recEn, U, m_meshObj,
			     uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg,  uPlusHalfPos);

    // cell jacobian functor
    // currently, REGARDLESS of the reconstruction scheme,
    // we have only first-order Jacobian so we need a first-order reconstructor
    // for the Jacobian
    scalar_type uMinusHalfNegForJ;
    scalar_type uMinusHalfPosForJ;
    scalar_type uPlusHalfNegForJ;
    scalar_type uPlusHalfPosForJ;
    reconstruct_functor_t ReconstructorForJ(::pressiodemoapps::InviscidFluxReconstruction::FirstOrder,
				 U, m_meshObj,
				 uMinusHalfNegForJ, uMinusHalfPosForJ,
				 uPlusHalfNegForJ,  uPlusHalfPosForJ);

    using jac_fnct_t = impl::FirstOrderInnerCellJacobianFunctor<
      base_t::dimensionality, numDofPerCell, scalar_type, jacobian_type, scalar_type, MeshType>;
    jac_fnct_t CellJacobianFunctor(m_jacobian, m_meshObj, JLneg, JLpos, JRneg, JRpos);

    // loop
    const auto dxInv = m_meshObj.dxInv();
    const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
    {
      Reconstructor.template operator()<numDofPerCell>(smPt);

      // note that REGARDLESS of the reconstruction scheme,
      // we currently only have only first-order Jacobian so we need
      // to run the reconstructor for the Jacobian
      // which will ensure that uMinusNegForJ, etc have the right values
      ReconstructorForJ.template operator()<numDofPerCell>(smPt);

      linAdvRusanovFlux(FL, uMinusHalfNeg, uMinusHalfPos);
      linAdvRusanovFlux(FR, uPlusHalfNeg,  uPlusHalfPos);
      linAdvRusanovFluxJacobian(JLneg, JLpos, uMinusHalfNegForJ, uMinusHalfPosForJ);
      linAdvRusanovFluxJacobian(JRneg, JRpos, uPlusHalfNegForJ,  uPlusHalfPosForJ);

      const auto vIndex = smPt*numDofPerCell;
      V(vIndex) = dxInv*(FL - FR);

      CellJacobianFunctor(smPt);
    }
  }

private:
  using base_t::m_meshObj;
  using base_t::m_recEn;
  using base_t::m_numDofStencilMesh;
  using base_t::m_numDofSampleMesh;
  using base_t::m_stencilVals;

  mutable jacobian_type m_jacobian = {};
  std::size_t m_jacNonZerosCount = {};
  bool m_onlyComputeVelocity = false;
};
#endif

}}//end namespace
#endif
