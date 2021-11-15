
#ifndef PRESSIODEMOAPPS_SWE2D_IMPL_HPP_
#define PRESSIODEMOAPPS_SWE2D_IMPL_HPP_

#include "swe_fluxes.hpp"
#include "swe_flux_jacobian.hpp"
#include "swe_2d_initial_condition.hpp"
#include "swe_2d_ghost_filler_inviscid_wall.hpp"
#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "cell_jacobian_first_order.hpp"
#include "Eigen/Sparse"

namespace pressiodemoapps{ namespace implswe{

template<class MeshType>
class EigenSwe2dApp
{

public:
  using index_t	      = typename MeshType::index_t;
  using scalar_type   = typename MeshType::scalar_t;
  using state_type    = Eigen::Matrix<scalar_type,-1,1>;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, index_t>;

  static constexpr int	   dimensionality{2};
  static constexpr index_t numDofPerCell{3};

private:
  using ghost_container_type   = Eigen::Matrix<scalar_type,-1,-1, Eigen::RowMajor>;
  using stencil_container_type = Eigen::Matrix<scalar_type,-1,1>;
  using flux_type	       = Eigen::Matrix<scalar_type,-1,1>;
  using edge_rec_type	       = Eigen::Matrix<scalar_type,-1,1>;
  using flux_jac_type = Eigen::Matrix<scalar_type, numDofPerCell, numDofPerCell>;

public:
  EigenSwe2dApp(const MeshType & meshObj,
		::pressiodemoapps::Swe2d probEn,
		::pressiodemoapps::InviscidFluxReconstruction recEn,
		::pressiodemoapps::InviscidFluxScheme fluxEnum,
		int icIdentifier)
    : m_icIdentifier(icIdentifier)
    ,m_probEn(probEn)
    ,m_recEn(recEn)
    ,m_fluxEn(fluxEnum)
    ,m_meshObj(meshObj)
  {
    m_numDofStencilMesh = m_meshObj.stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * numDofPerCell;
    allocateStencilValuesContainer();
    allocateGhosts();
  }

  scalar_type gravity()         const{ return m_gravity; }
  scalar_type coriolis()        const{ return m_coriolis; }

  state_type initialCondition() const{
    return initialConditionImpl();
  }

protected:
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

  void initializeJacobian(jacobian_type & J)
  {
    J.resize(m_numDofSampleMesh, m_numDofStencilMesh);
    initializeJacobianFirstOrder(J);
    if (!J.isCompressed()){
      J.makeCompressed();
    }
  }

  void initializeJacobianFirstOrder(jacobian_type & J)
  {
    using Tr = Eigen::Triplet<scalar_type>;
    std::vector<Tr> trList;

    const scalar_type val0 = 0.;
    const auto & graph = m_meshObj.graph();
    for (int cell=0; cell<m_meshObj.sampleMeshSize(); ++cell)
      {
	const auto jacRow = cell*numDofPerCell;
	const auto ci0  = graph(cell, 0)*numDofPerCell;

	for (int k=0; k<numDofPerCell; ++k){
	  for (int j=0; j<numDofPerCell; ++j){
	    trList.push_back( Tr(jacRow+k, ci0+j, val0) );
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
	      trList.push_back( Tr(jacRow+k, ci+j, val0) );
	    }
	  }
	}

	if (F0 != -1){
	  const auto ci = F0*numDofPerCell;
	  for (int k=0; k<numDofPerCell; ++k){
	    for (int j=0; j<numDofPerCell; ++j){
	      trList.push_back( Tr(jacRow+k, ci+j, val0) );
	    }
	  }
	}

	if (R0 != -1){
	  const auto ci = R0*numDofPerCell;
	  for (int k=0; k<numDofPerCell; ++k){
	    for (int j=0; j<numDofPerCell; ++j){
	      trList.push_back( Tr(jacRow+k, ci+j, val0) );
	    }
	  }
	}

	if (B0 != -1){
	  const auto ci = B0*numDofPerCell;
	  for (int k=0; k<numDofPerCell; ++k){
	    for (int j=0; j<numDofPerCell; ++j){
	      trList.push_back( Tr(jacRow+k, ci+j, val0) );
	    }
	  }
	}
      }

    J.setFromTriplets(trList.begin(), trList.end());
  }

  template<class U_t>
  void fillGhosts(const U_t & U, const scalar_type currentTime) const
  {
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    if (m_probEn == ::pressiodemoapps::Swe2d::SlipWall)
    {
      using ghost_filler_t  = ::pressiodemoapps::implswe::InviscidWallFiller<
	U_t, MeshType, ghost_container_type>;
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


  // note that here we MUST use a template because when doing
  // bindings, this gets deduced to be a Ref
  template<class U_t, class V_t>
  void velocityAndOptionalJacobian(const U_t & U,
				   const scalar_type currentTime,
				   V_t & V,
				   jacobian_type * J) const
  {
    flux_type FL(numDofPerCell);
    flux_type FR(numDofPerCell);
    flux_type FU(numDofPerCell);
    flux_type FD(numDofPerCell);
    edge_rec_type uMinusHalfNeg(numDofPerCell);
    edge_rec_type uMinusHalfPos(numDofPerCell);
    edge_rec_type uPlusHalfNeg (numDofPerCell);
    edge_rec_type uPlusHalfPos (numDofPerCell);

    fillGhosts(U, currentTime);

    if (J){
      ::pressiodemoapps::set_zero(*J);
    }

    velocityCellsNearBdImpl(U, currentTime, V, J,
			    FL, FR, FD, FU,
			    uMinusHalfNeg, uMinusHalfPos,
			    uPlusHalfNeg,  uPlusHalfPos);

    velocityInnerCellsImpl(U, currentTime, V, J,
			   FL, FR, FD, FU,
			   uMinusHalfNeg, uMinusHalfPos,
			   uPlusHalfNeg,  uPlusHalfPos);
  }

  void allocateStencilValuesContainer()
  {
    // the stencil size needed is determined by the desired reconstruction
    // kind NOT from the mesh. THis is important because for example
    // the mesh can have a wider connectivity that what is needed.
    const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_recEn);
    ::pressiodemoapps::resize(m_stencilVals, numDofPerCell*stencilSizeNeeded);
  }

  template<class U_t, class V_t>
  void velocityCellsNearBdImpl(const U_t & U,
			       const scalar_type currentTime,
			       V_t & V,
			       jacobian_type * J,
			       flux_type & FL,
			       flux_type & FR,
			       flux_type & FD,
			       flux_type & FU,
			       edge_rec_type & uMinusHalfNeg,
			       edge_rec_type & uMinusHalfPos,
			       edge_rec_type & uPlusHalfNeg,
			       edge_rec_type & uPlusHalfPos) const
  {

    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    const auto dxInv   = m_meshObj.dxInv();
    const auto dyInv   = m_meshObj.dyInv();

    using stencil_filler_t = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_container_type, U_t, MeshType, ghost_container_type>;
    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromStencilThreeDofPerCell
      <edge_rec_type, stencil_container_type>;

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

    // ----------------------------------------
    // functors needed to compute cell jacobian
    // we have only first-order Jacobian for now, so we need
    // dedicate functors because we cannot use those above.
    // Once we have cell Jacobians of various order, we can change this.
    // ----------------------------------------
    flux_jac_type JLneg, JLpos;
    flux_jac_type JRneg, JRpos;
    flux_jac_type JFneg, JFpos;
    flux_jac_type JBneg, JBpos;

    const auto cellJacOrdEn = ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder;

    stencil_container_type stencilValsForJ = {};
    const auto stencilSizeForJ = reconstructionTypeToStencilSize(cellJacOrdEn);
    ::pressiodemoapps::resize(stencilValsForJ, numDofPerCell*stencilSizeForJ);
    stencil_filler_t FillStencilValuesXFunctorForJ(stencilSizeForJ, U, m_meshObj,
						   m_ghostLeft, m_ghostRight, stencilValsForJ,
						   xAxis);
    stencil_filler_t FillStencilValuesYFunctorForJ(stencilSizeForJ, U, m_meshObj,
						   m_ghostBack, m_ghostFront, stencilValsForJ,
						   yAxis);

    edge_rec_type uMinusHalfNegForJ(numDofPerCell);
    edge_rec_type uMinusHalfPosForJ(numDofPerCell);
    edge_rec_type uPlusHalfNegForJ(numDofPerCell);
    edge_rec_type uPlusHalfPosForJ(numDofPerCell);
    rec_fnct_t FaceValuesReconstructFunctorForJ(cellJacOrdEn, stencilValsForJ,
						uMinusHalfNegForJ, uMinusHalfPosForJ,
						uPlusHalfNegForJ,  uPlusHalfPosForJ);

    using jac_fnct_t = ::pressiodemoapps::impl::FirstOrderBdCellJacobianFunctor<
      dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, MeshType>;

    jac_fnct_t CellJacobianFunctorX(*J, m_meshObj, JLneg, JLpos, JRneg, JRpos, xAxis);
    jac_fnct_t CellJacobianFunctorY(*J, m_meshObj, JBneg, JBpos, JFneg, JFpos, yAxis);

    std::array<scalar_type, numDofPerCell> bcCellJacFactorsX;
    std::array<scalar_type, numDofPerCell> bcCellJacFactorsY;
    bcCellJacFactorsX.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsY.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsX[1] = static_cast<scalar_type>(-1);
    bcCellJacFactorsY[2] = static_cast<scalar_type>(-1);

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

      if (J){
	FillStencilValuesXFunctorForJ(smPt, it);
	FaceValuesReconstructFunctorForJ();
      }

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  sweRusanovFluxThreeDof(FL, uMinusHalfNeg, uMinusHalfPos, normalX_, m_gravity);
	  sweRusanovFluxThreeDof(FR, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gravity);

	  if (J){
	    sweRusanovFluxJacobianThreeDof(JLneg, JLpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
					   normalX_, m_gravity);
	    sweRusanovFluxJacobianThreeDof(JRneg, JRpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
					   normalX_, m_gravity);
	  }

	  break;
	}

      if (J){
	CellJacobianFunctorX(smPt, bcCellJacFactorsX, 1);
      }

      // *** Y ***
      StencilFillerY(smPt, it);
      Reconstructor();

      if (J){
	FillStencilValuesYFunctorForJ(smPt, it);
	FaceValuesReconstructFunctorForJ();
      }

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  sweRusanovFluxThreeDof(FD, uMinusHalfNeg, uMinusHalfPos, normalY_, m_gravity);
	  sweRusanovFluxThreeDof(FU, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gravity);

	  if (J){
	    sweRusanovFluxJacobianThreeDof(JBneg, JBpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
					   normalY_, m_gravity);
	    sweRusanovFluxJacobianThreeDof(JFneg, JFpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
					   normalY_, m_gravity);
	  }

	  break;
	}

      if (J){
	CellJacobianFunctorY(smPt, bcCellJacFactorsY, 1);
      }

      const auto vIndex = smPt*numDofPerCell;
      const auto uIndex = graph(smPt, 0)*numDofPerCell;
      V(vIndex)   = dxInv*(FL(0) - FR(0)) + dyInv*(FD(0) - FU(0));
      V(vIndex+1) = dxInv*(FL(1) - FR(1)) + dyInv*(FD(1) - FU(1)) - m_coriolis*U(uIndex+2)/U(uIndex);
      V(vIndex+2) = dxInv*(FL(2) - FR(2)) + dyInv*(FD(2) - FU(2)) + m_coriolis*U(uIndex+1)/U(uIndex);

      if (J){
	const index_t col_i = graph(smPt, 0)*numDofPerCell;
	J->coeffRef(vIndex+1, col_i)   +=  m_coriolis*U(uIndex+2)/(U(uIndex)*U(uIndex) );
	J->coeffRef(vIndex+1, col_i+2) += -m_coriolis/U(uIndex);
	J->coeffRef(vIndex+2, col_i+1) +=  m_coriolis/U(uIndex);
	J->coeffRef(vIndex+2, col_i)   +=  -m_coriolis*U(uIndex+1)/(U(uIndex)*U(uIndex) );
      }
    }
  }

  template<class U_t, class V_t>
  void velocityInnerCellsImpl(const U_t & U,
			      const scalar_type currentTime,
			      V_t & V,
			      jacobian_type * J,
			      flux_type & FL,
			      flux_type & FR,
			      flux_type & FD,
			      flux_type & FU,
			      edge_rec_type & uMinusHalfNeg,
			      edge_rec_type & uMinusHalfPos,
			      edge_rec_type & uPlusHalfNeg,
			      edge_rec_type & uPlusHalfPos) const
  {

    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    const auto dxInv   = m_meshObj.dxInv();
    const auto dyInv   = m_meshObj.dyInv();

    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromState<
      dimensionality, edge_rec_type, U_t, MeshType>;

    using jac_fnct_t = ::pressiodemoapps::impl::FirstOrderInnerCellJacobianFunctor<
      dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, MeshType>;

    // reconstruct functor for face fluxes
    // here we need to use whatever order (m_recEn) user decides
    rec_fnct_t ReconstructorX(xAxis, m_recEn, U, m_meshObj,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);

    rec_fnct_t ReconstructorY(yAxis, m_recEn, U, m_meshObj,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);

    // cell jacobian functor
    // currently, REGARDLESS of the reconstruction scheme,
    // we have only first-order Jacobian so we need a first-order reconstructor for the Jacobian
    flux_jac_type JLneg, JLpos;
    flux_jac_type JRneg, JRpos;
    flux_jac_type JFneg, JFpos;
    flux_jac_type JBneg, JBpos;

    const auto firstOrder = ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder;
    edge_rec_type uMinusHalfNegForJ(numDofPerCell);
    edge_rec_type uMinusHalfPosForJ(numDofPerCell);
    edge_rec_type uPlusHalfNegForJ(numDofPerCell);
    edge_rec_type uPlusHalfPosForJ(numDofPerCell);
    rec_fnct_t ReconstructorXForJ(xAxis, firstOrder, U, m_meshObj,
				  uMinusHalfNegForJ, uMinusHalfPosForJ,
				  uPlusHalfNegForJ,  uPlusHalfPosForJ);

    rec_fnct_t ReconstructorYForJ(yAxis, firstOrder, U, m_meshObj,
				  uMinusHalfNegForJ, uMinusHalfPosForJ,
				  uPlusHalfNegForJ,  uPlusHalfPosForJ);

    jac_fnct_t CellJacobianFunctorX(*J, m_meshObj, JLneg, JLpos, JRneg, JRpos, xAxis);
    jac_fnct_t CellJacobianFunctorY(*J, m_meshObj, JBneg, JBpos, JFneg, JFpos, yAxis);

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

      if (J){
	// note that REGARDLESS of the reconstruction scheme,
	// we currently only have only first-order Jacobian so we need
	// to run the reconstructor for the Jacobian
	// which will ensure that uMinusNegForJ, etc have the right values
	ReconstructorXForJ.template operator()<numDofPerCell>(smPt);
      }

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  sweRusanovFluxThreeDof(FL, uMinusHalfNeg, uMinusHalfPos, normalX_, m_gravity);
	  sweRusanovFluxThreeDof(FR, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gravity);

	  if (J){
	    sweRusanovFluxJacobianThreeDof(JLneg, JLpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
					   normalX_, m_gravity);
	    sweRusanovFluxJacobianThreeDof(JRneg, JRpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
					 normalX_, m_gravity);
	  }
	  break;
	}

      if (J){
	CellJacobianFunctorX(smPt);
      }

      // Y
      ReconstructorY.template operator()<numDofPerCell>(smPt);

      if (J){
	// note that REGARDLESS of the reconstruction scheme,
	// we currently only have only first-order Jacobian so we need
	// to run the reconstructor for the Jacobian
	// which will ensure that uMinusNegForJ, etc have the right values
	ReconstructorYForJ.template operator()<numDofPerCell>(smPt);
      }

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  sweRusanovFluxThreeDof(FD, uMinusHalfNeg, uMinusHalfPos, normalY_, m_gravity);
	  sweRusanovFluxThreeDof(FU, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gravity);

	  if (J){
	    sweRusanovFluxJacobianThreeDof(JBneg, JBpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
					   normalY_, m_gravity);
	    sweRusanovFluxJacobianThreeDof(JFneg, JFpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
					   normalY_, m_gravity);
	  }
	  break;
	}

      if (J){
	CellJacobianFunctorY(smPt);
      }

      V(vIndex)   = dxInv*(FL(0) - FR(0)) + dyInv*(FD(0) - FU(0));
      V(vIndex+1) = dxInv*(FL(1) - FR(1)) + dyInv*(FD(1) - FU(1)) - m_coriolis*U(uIndex+2)/U(uIndex);
      V(vIndex+2) = dxInv*(FL(2) - FR(2)) + dyInv*(FD(2) - FU(2)) + m_coriolis*U(uIndex+1)/U(uIndex);

      if (J){
	const index_t col_i = graph(smPt, 0)*numDofPerCell;
	J->coeffRef(vIndex+1, col_i)   +=  m_coriolis*U(uIndex+2)/(U(uIndex)*U(uIndex) );
	J->coeffRef(vIndex+1, col_i+2) += -m_coriolis/U(uIndex);
	J->coeffRef(vIndex+2, col_i+1) +=  m_coriolis/U(uIndex);
	J->coeffRef(vIndex+2, col_i)   +=  -m_coriolis*U(uIndex+1)/(U(uIndex)*U(uIndex) );
      }
    }
  }


private:
  void allocateGhosts()
  {
    const auto stencilSize    = reconstructionTypeToStencilSize(m_recEn);
    const auto numGhostValues = numDofPerCell*((stencilSize-1)/2);

    const index_t s1 = m_meshObj.numCellsBd();
    ::pressiodemoapps::resize(m_ghostLeft, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostFront,s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight,s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostBack, s1, numGhostValues);
  }

protected:
  const int m_icIdentifier;
  const scalar_type m_gravity  = static_cast<scalar_type>(9.8);
  const scalar_type m_coriolis = static_cast<scalar_type>(-3.0);
  const scalar_type m_initialPulseMagnitude = static_cast<scalar_type>(0.125);

  ::pressiodemoapps::Swe2d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  ::pressiodemoapps::InviscidFluxScheme m_fluxEn;

  const MeshType & m_meshObj;
  mutable stencil_container_type m_stencilVals;

  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  mutable ghost_container_type m_ghostLeft;
  mutable ghost_container_type m_ghostFront;
  mutable ghost_container_type m_ghostRight;
  mutable ghost_container_type m_ghostBack;

  const std::array<scalar_type, 2> normalX_{1, 0};
  const std::array<scalar_type, 2> normalY_{0, 1};
};

}}//end namespace
#endif
