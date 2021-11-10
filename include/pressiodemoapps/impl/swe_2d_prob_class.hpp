
#ifndef PRESSIODEMOAPPS_SWE2D_IMPL_HPP_
#define PRESSIODEMOAPPS_SWE2D_IMPL_HPP_

// note that the code below is pretty ugly for now, but it works.
// this is inside impl namespace for a reason, and will need to be
// improved later on but we have a starting point.

#include "swe_fluxes.hpp"
#include "swe_2d_initial_condition.hpp"
#include "swe_2d_ghost_filler_inviscid_wall.hpp"
#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "cell_jacobian_first_order.hpp"

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
#include "Eigen/Sparse"
#endif

namespace pressiodemoapps{ namespace implswe{

template<
  class scalar_t,
  class mesh_t,
  class state_t,
  class velo_t,
  class ghost_t,
  class jacobian_t = void
  >
class Swe2dAppT
{

public:
  using index_t		 = typename mesh_t::index_t;
  using scalar_type	 = scalar_t;
  using state_type	 = state_t;
  using velocity_type	 = state_t;
  using ghost_type	 = ghost_t;
  using stencil_values_t = state_type;
  using flux_t		 = state_type;
  using edge_rec_t	 = state_type;

  static constexpr int dimensionality{2};
  static constexpr index_t numDofPerCell{3};

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
  using jacobian_type = jacobian_t;
  using flux_jac_type = Eigen::Matrix<scalar_type, numDofPerCell, numDofPerCell>;
#endif

public:
#if not defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  Swe2dAppT(const mesh_t & meshObj,
	      ::pressiodemoapps::Swe2d probEn,
	      ::pressiodemoapps::InviscidFluxReconstruction recEn,
	      ::pressiodemoapps::InviscidFluxScheme fluxEnum,
        int icIdentifier)
    : m_icIdentifier(icIdentifier),
      m_probEn(probEn),
      m_recEn(recEn),
      m_fluxEn(fluxEnum),
      m_meshObj(meshObj)
  {

    computeDofs();
    allocateStencilValuesContainer();
    allocateGhosts();
    initializeJacobian();
  }

#else
  // note that when doing bindings, I need to first construct
  // ghost with {1,1} just so that the numpy array picks up they
  // are 2dim array otherwise it thinks they are 1d array.
  // The right allocation for these is then done inside allocateGhosts.
  Swe2dAppT(const mesh_t & meshObj,
	      ::pressiodemoapps::Swe2d probEn,
	      ::pressiodemoapps::InviscidFluxReconstruction recEn,
	      ::pressiodemoapps::InviscidFluxScheme fluxEnum,
        int icIdentifier)
    : m_icIdentifier(icIdentifier),
      m_probEn(probEn),
      m_recEn(recEn),
      m_fluxEn(fluxEnum),
      m_meshObj(meshObj),
      m_stencilVals(1),
      m_ghostLeft({1,1}),
      m_ghostFront({1,1}),
      m_ghostRight({1,1}),
      m_ghostBack({1,1})
  {
    computeDofs();
    allocateStencilValuesContainer();
    allocateGhosts();
  }
#endif

  scalar_type gravity()         const{ return m_gravity; }
  scalar_type coriolis()        const{ return m_coriolis; }
  index_t totalDofSampleMesh()  const{ return m_numDofSampleMesh; }
  index_t totalDofStencilMesh() const{ return m_numDofStencilMesh; }

  state_type initialCondition() const{
    return initialConditionImpl();
  }

  // the Jacobian is by default fused with the velocity,
  // this method allows one to disable the jacobian
  // so only velocity is computed
  void disableJacobian() {
    m_onlyComputeVelocity = true;
  }

  velocity_type createVelocity() const{
    velocity_type V(m_numDofSampleMesh);
    return V;
  }

  jacobian_type createJacobian() const {
    return m_jacobian;
  }

  void velocity(const state_type & state,
		const scalar_type currentTime,
		velocity_type & V) const
  {
    flux_t FL(numDofPerCell);
    flux_t FR(numDofPerCell);
    flux_t FU(numDofPerCell);
    flux_t FD(numDofPerCell);
    edge_rec_t uMinusHalfNeg(numDofPerCell);
    edge_rec_t uMinusHalfPos(numDofPerCell);
    edge_rec_t uPlusHalfNeg (numDofPerCell);
    edge_rec_t uPlusHalfPos (numDofPerCell);

    fillGhosts(state, currentTime);
    velocityCellsNearBdImpl(state, currentTime, V, FL, FR, FD, FU,
			    uMinusHalfNeg, uMinusHalfPos,
			    uPlusHalfNeg,  uPlusHalfPos);
    velocityInnerCellsImpl(state, currentTime, V, FL, FR, FD, FU,
			   uMinusHalfNeg, uMinusHalfPos,
			   uPlusHalfNeg,  uPlusHalfPos);
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
  void computeDofs(){
    // calculate total num of dofs on sample and stencil mesh
    m_numDofStencilMesh = m_meshObj.stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * numDofPerCell;
  }

  void allocateStencilValuesContainer()
  {
    // the stencil size needed is determined by the desired reconstruction
    // kind NOT from the mesh. THis is important because for example
    // the mesh can have a wider connectivity that what is needed.
    const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_recEn);
    ::pressiodemoapps::resize(m_stencilVals, numDofPerCell*stencilSizeNeeded);
  }

  state_type initialConditionImpl() const
  {
    state_type initialState(m_numDofStencilMesh);

    switch(m_probEn)
      {
      case ::pressiodemoapps::Swe2d::SlipWall:{
	if( m_icIdentifier == 1){
	  GaussianPulse(initialState, m_meshObj, m_initialPulseMagnitude);
        }
        else{
	  throw std::runtime_error("invalid IC");
        }
      }
      };

    return initialState;
  }

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
  void initializeJacobian()
  {
    m_jacobian.resize(m_numDofSampleMesh, m_numDofStencilMesh);
    initializeJacobianFirstOrder();

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
    using Tr = Eigen::Triplet<scalar_type>;
    std::vector<Tr> trList;

    const scalar_type val0 = 0.;
    const auto & graph = m_meshObj.graph();
    for (int cell=0; cell<m_meshObj.sampleMeshSize(); ++cell)
      {
	const auto jacRowOfCurrentCellDensity = cell*numDofPerCell;
	const auto ci0  = graph(cell, 0)*numDofPerCell;

	for (int k=0; k<numDofPerCell; ++k){
	  for (int j=0; j<numDofPerCell; ++j){
	    trList.push_back( Tr(jacRowOfCurrentCellDensity+k, ci0+j, val0) );
	  }
	}

	const auto L0 = graph(cell, 1);
	const auto F0 = graph(cell, 2);
	const auto R0 = graph(cell, 3);
	const auto B0 = graph(cell, 4);
	if (L0 != -1){
	  const auto ci = L0*numDofPerCell;
	  for (int k=0; k<numDofPerCell; ++k){
	    for (int j=0; j<numDofPerCell; ++j){
	      trList.push_back( Tr(jacRowOfCurrentCellDensity+k, ci+j, val0) );
	    }
	  }
	}

	if (F0 != -1){
	  const auto ci = F0*numDofPerCell;
	  for (int k=0; k<numDofPerCell; ++k){
	    for (int j=0; j<numDofPerCell; ++j){
	      trList.push_back( Tr(jacRowOfCurrentCellDensity+k, ci+j, val0) );
	    }
	  }
	}

	if (R0 != -1){
	  const auto ci = R0*numDofPerCell;
	  for (int k=0; k<numDofPerCell; ++k){
	    for (int j=0; j<numDofPerCell; ++j){
	      trList.push_back( Tr(jacRowOfCurrentCellDensity+k, ci+j, val0) );
	    }
	  }
	}

	if (B0 != -1){
	  const auto ci = B0*numDofPerCell;
	  for (int k=0; k<numDofPerCell; ++k){
	    for (int j=0; j<numDofPerCell; ++j){
	      trList.push_back( Tr(jacRowOfCurrentCellDensity+k, ci+j, val0) );
	    }
	  }
	}
      }

    m_jacobian.setFromTriplets(trList.begin(), trList.end());
  }

#endif

  void velocityCellsNearBdImpl(const state_type & U,
			       const scalar_type currentTime,
			       velocity_type & V,
			       flux_t & FL,
			       flux_t & FR,
			       flux_t & FD,
			       flux_t & FU,
			       edge_rec_t & uMinusHalfNeg,
			       edge_rec_t & uMinusHalfPos,
			       edge_rec_t & uPlusHalfNeg,
			       edge_rec_t & uPlusHalfPos) const
  {

    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    const auto dxInv   = m_meshObj.dxInv();
    const auto dyInv   = m_meshObj.dyInv();

    flux_jac_type JLneg, JLpos;
    flux_jac_type JRneg, JRpos;
    flux_jac_type JFneg, JFpos;
    flux_jac_type JBneg, JBpos;

    using stencil_filler_t = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_values_t, state_type, mesh_t, ghost_t>;
    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromStencilThreeDofPerCell
      <edge_rec_t, stencil_values_t>;

    using jac_fnct_t = ::pressiodemoapps::impl::FirstOrderBdCellJacobianFunctor<
      dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, mesh_t>;

    // stencil filler and reconstructor for face fluxes
    // here we need to use whatever order (m_recEn) user decides
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    stencil_filler_t StencilFillerX(stencilSize, U, m_meshObj,
				    m_ghostLeft, m_ghostRight,
				    m_stencilVals, xAxis);

    stencil_filler_t StencilFillerY(stencilSize, U, m_meshObj,
				    m_ghostBack, m_ghostFront,
				    m_stencilVals, yAxis);

    rec_fnct_t Reconstructor(m_recEn, m_stencilVals,
			     uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg,  uPlusHalfPos);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
    // ----------------------------------------
    // functors needed to compute cell jacobian
    // we have only first-order Jacobian for now, so we need
    // dedicate functors because we cannot use those above.
    // Once we have cell Jacobians of various order, we can change this.
    // ----------------------------------------
    const auto cellJacOrdEn = ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder;

    stencil_values_t stencilValsForJ = {};
    const auto stencilSizeForJ = reconstructionTypeToStencilSize(cellJacOrdEn);
    ::pressiodemoapps::resize(stencilValsForJ, numDofPerCell*stencilSizeForJ);
    stencil_filler_t FillStencilValuesXFunctorForJ(stencilSizeForJ, U, m_meshObj,
						   m_ghostLeft, m_ghostRight, stencilValsForJ, xAxis);
    stencil_filler_t FillStencilValuesYFunctorForJ(stencilSizeForJ, U, m_meshObj,
						   m_ghostBack, m_ghostFront, stencilValsForJ, yAxis);

    edge_rec_t uMinusHalfNegForJ(numDofPerCell);
    edge_rec_t uMinusHalfPosForJ(numDofPerCell);
    edge_rec_t uPlusHalfNegForJ(numDofPerCell);
    edge_rec_t uPlusHalfPosForJ(numDofPerCell);
    rec_fnct_t FaceValuesReconstructFunctorForJ(cellJacOrdEn, stencilValsForJ,
						uMinusHalfNegForJ, uMinusHalfPosForJ,
						uPlusHalfNegForJ,  uPlusHalfPosForJ);

    jac_fnct_t CellJacobianFunctorX(m_jacobian, m_meshObj, JLneg, JLpos, JRneg, JRpos, xAxis);
    jac_fnct_t CellJacobianFunctorY(m_jacobian, m_meshObj, JBneg, JBpos, JFneg, JFpos, yAxis);
#endif

    // ----------------------------------------
    // loop over cells
    // ----------------------------------------
    const auto & graph = m_meshObj.graph();
    const auto & rows = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<rows.size(); ++it)
    {
      const auto smPt = rows[it];

      // X
      StencilFillerX(smPt, it);
      Reconstructor();

      if (!m_onlyComputeVelocity){
	FillStencilValuesXFunctorForJ(smPt, it);
	FaceValuesReconstructFunctorForJ();
      }

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  sweRusanovFluxThreeDof(FL, uMinusHalfNeg, uMinusHalfPos, normalX_, m_gravity);
	  sweRusanovFluxThreeDof(FR, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gravity);

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
// 	  if (!m_onlyComputeVelocity){
// 	    eeRusanovFluxJacobianFourDof(JLneg, JLpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
// 					 normalX_, m_gravity);
// 	    eeRusanovFluxJacobianFourDof(JRneg, JRpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
// 					 normalX_, m_gravity);
// 	  }
// #endif

	  break;
	}

      // Y
      StencilFillerY(smPt, it);
      Reconstructor();

      if (!m_onlyComputeVelocity){
	FillStencilValuesYFunctorForJ(smPt, it);
	FaceValuesReconstructFunctorForJ();
      }

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  sweRusanovFluxThreeDof(FD, uMinusHalfNeg, uMinusHalfPos, normalY_, m_gravity);
	  sweRusanovFluxThreeDof(FU, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gravity);

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
// 	  if (!m_onlyComputeVelocity){
// 	    eeRusanovFluxJacobianFourDof(JBneg, JBpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
// 					 normalY_, m_gravity);
// 	    eeRusanovFluxJacobianFourDof(JFneg, JFpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
// 					 normalY_, m_gravity);
// 	  }
// #endif

	  break;
	}

      // if (!m_onlyComputeVelocity){
      // 	CellJacobianFunctorY(smPt, bdType);
      // }

      const auto vIndex = smPt*numDofPerCell;
      const auto uIndex = graph(smPt, 0)*numDofPerCell;

      V(vIndex)   = dxInv*(FL(0) - FR(0)) + dyInv*(FD(0) - FU(0));
      V(vIndex+1) = dxInv*(FL(1) - FR(1)) + dyInv*(FD(1) - FU(1)) - m_coriolis*U(uIndex+2)/U(uIndex);
      V(vIndex+2) = dxInv*(FL(2) - FR(2)) + dyInv*(FD(2) - FU(2)) + m_coriolis*U(uIndex+1)/U(uIndex);
    }
  }

  void velocityInnerCellsImpl(const state_type & U,
			      const scalar_type currentTime,
			      velocity_type & V,
			      flux_t & FL,
			      flux_t & FR,
			      flux_t & FD,
			      flux_t & FU,
			      edge_rec_t & uMinusHalfNeg,
			      edge_rec_t & uMinusHalfPos,
			      edge_rec_t & uPlusHalfNeg,
			      edge_rec_t & uPlusHalfPos) const
  {

    constexpr int xAxis = 1;
    constexpr int yAxis = 2;

    const auto dxInv   = m_meshObj.dxInv();
    const auto dyInv   = m_meshObj.dyInv();

    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromState<
      dimensionality, edge_rec_t, state_type, mesh_t>;

    rec_fnct_t ReconstructorX(xAxis, m_recEn, U, m_meshObj,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);

    rec_fnct_t ReconstructorY(yAxis, m_recEn, U, m_meshObj,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);

    // deal with cells away from boundaries
    // const auto & graph = m_meshObj.graph();
    const auto & graph = m_meshObj.graph();
    const auto & rowsIn = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (int it=0; it<rowsIn.size(); ++it)
    {
      const auto smPt = rowsIn[it];
      const auto vIndex = smPt*numDofPerCell;
      const auto uIndex = graph(smPt, 0)*numDofPerCell;

      // X
      ReconstructorX.template operator()<numDofPerCell>(smPt);
      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  sweRusanovFluxThreeDof(FL, uMinusHalfNeg, uMinusHalfPos, normalX_, m_gravity);
	  sweRusanovFluxThreeDof(FR, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gravity);
	  break;
	}

      // Y
      ReconstructorY.template operator()<numDofPerCell>(smPt);
      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  sweRusanovFluxThreeDof(FD, uMinusHalfNeg, uMinusHalfPos, normalY_, m_gravity);
	  sweRusanovFluxThreeDof(FU, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gravity);
	  break;
	}


      V(vIndex)   = dxInv*(FL(0) - FR(0)) + dyInv*(FD(0) - FU(0));
      V(vIndex+1) = dxInv*(FL(1) - FR(1)) + dyInv*(FD(1) - FU(1)) - m_coriolis*U(uIndex+2)/U(uIndex);
      V(vIndex+2) = dxInv*(FL(2) - FR(2)) + dyInv*(FD(2) - FU(2)) + m_coriolis*U(uIndex+1)/U(uIndex);
    }

  }

  void fillGhosts(const state_type & U, const scalar_type currentTime) const
  {
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    if (m_probEn == ::pressiodemoapps::Swe2d::SlipWall)
    {
      using ghost_filler_t  = ::pressiodemoapps::implswe::InviscidWallFiller<
	state_type, mesh_t, ghost_t>;
      ghost_filler_t ghF(stencilSize, U, m_meshObj,
			 m_ghostLeft, m_ghostFront,
			 m_ghostRight, m_ghostBack);

      const auto & rowsBd = m_meshObj.graphRowsOfCellsNearBd();
      if (stencilSize==3){
	for (int it=0; it<rowsBd.size(); ++it){
	  ghF.template operator()<3>(rowsBd[it], it);
	}
      }
      else if(stencilSize==5){
	for (int it=0; it<rowsBd.size(); ++it){
	  ghF.template operator()<5>(rowsBd[it], it);
	}
      }
      else{
	for (int it=0; it<rowsBd.size(); ++it){
	  ghF.template operator()<7>(rowsBd[it], it);
	}
      }
    }

    else{
      // no op
    }

  }

private:
  void allocateGhosts()
  {
    const auto stencilSize    = reconstructionTypeToStencilSize(m_recEn);
    const auto numGhostValues = numDofPerCell*((stencilSize-1)/2);

    const index_t s1 = m_meshObj.numCellsBd();
    ::pressiodemoapps::resize(m_ghostLeft,   s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostFront,    s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight,  s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostBack, s1, numGhostValues);
  }

private:
  const int m_icIdentifier;
  const scalar_type m_gravity  = static_cast<scalar_type>(9.8);
  const scalar_type m_coriolis = static_cast<scalar_type>(-3.0);
  const scalar_type m_initialPulseMagnitude = static_cast<scalar_type>(0.125);

  ::pressiodemoapps::Swe2d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  ::pressiodemoapps::InviscidFluxScheme m_fluxEn;

  const mesh_t & m_meshObj;
  mutable stencil_values_t m_stencilVals;

  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  mutable ghost_t m_ghostLeft;
  mutable ghost_t m_ghostFront;
  mutable ghost_t m_ghostRight;
  mutable ghost_t m_ghostBack;

  const std::array<scalar_type, 2> normalX_{1, 0};
  const std::array<scalar_type, 2> normalY_{0, 1};

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
  mutable jacobian_type m_jacobian = {};
  std::size_t m_jacNonZerosCount = {};
  bool m_onlyComputeVelocity = false;
#endif
};

}}//end namespace
#endif
