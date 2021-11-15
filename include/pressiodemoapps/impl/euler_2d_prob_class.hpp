
#ifndef PRESSIODEMOAPPS_EULER2D_APP_HPP_
#define PRESSIODEMOAPPS_EULER2D_APP_HPP_

// note that the code below is pretty ugly for now, but it works.
// this is inside impl namespace for a reason, and will need to be
// improved later on but we have a starting point.

#include "euler_fluxes.hpp"
#include "euler_flux_jacobian.hpp"
#include "euler_rankine_hugoniot.hpp"
#include "euler_2d_initial_condition.hpp"
#include "functor_ghost_fill_neumann.hpp"
#include "euler_2d_ghost_filler_sedov2d_sym.hpp"
#include "euler_2d_ghost_filler_normal_shock.hpp"
#include "euler_2d_ghost_filler_double_mach_reflection.hpp"
#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "cell_jacobian_first_order.hpp"
#include "Eigen/Sparse"

namespace pressiodemoapps{ namespace ee{ namespace impl{

template<class MeshType>
class EigenEuler2dApp
{

public:
  using index_t	      = typename MeshType::index_t;
  using scalar_type   = typename MeshType::scalar_t;
  using state_type    = Eigen::Matrix<scalar_type,-1,1>;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, index_t>;

  static constexpr int	   dimensionality{2};
  static constexpr index_t numDofPerCell{4};

private:
  using ghost_container_type   = Eigen::Matrix<scalar_type,-1,-1, Eigen::RowMajor>;
  using stencil_container_type = Eigen::Matrix<scalar_type,-1,1>;
  using flux_type	       = Eigen::Matrix<scalar_type,-1,1>;
  using edge_rec_type	       = Eigen::Matrix<scalar_type,-1,1>;
  using flux_jac_type = Eigen::Matrix<scalar_type, numDofPerCell, numDofPerCell>;

public:
  EigenEuler2dApp(const MeshType & meshObj,
		  ::pressiodemoapps::Euler2d probEnum,
		  ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		  ::pressiodemoapps::InviscidFluxScheme fluxEnum,
		  int icIdentifier)
    : m_meshObj(meshObj), m_probEn(probEnum), m_recEn(recEnum),
      m_fluxEn(fluxEnum), m_icIdentifier(icIdentifier)
  {
    // calculate total num of dofs on sample and stencil mesh
    m_numDofStencilMesh = m_meshObj.stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * numDofPerCell;

    allocateStencilValuesContainer();
    allocateGhosts();
  }

  scalar_type gamma() const{
    return m_gamma;
  }

  state_type initialCondition() const{
    return initialConditionImpl();
  }

protected:
  state_type initialConditionImpl() const
  {
    state_type initialState(m_numDofStencilMesh);

    switch(m_probEn)
      {
      case ::pressiodemoapps::Euler2d::PeriodicSmooth:{
	sin2dEulerIC(initialState, m_meshObj, m_gamma);
	return initialState;
      }

      case ::pressiodemoapps::Euler2d::KelvinHelmholtz:{
	KelvinHelmholtzIC(initialState, m_meshObj, m_gamma);
	return initialState;
      }

      case ::pressiodemoapps::Euler2d::SedovFull:{
	sedov2dIC(initialState, m_meshObj, m_gamma);
	return initialState;
      }

      case ::pressiodemoapps::Euler2d::SedovSymmetry:{
	sedov2dsymmetryIC(initialState, m_meshObj, m_gamma);
	return initialState;
      }

      case ::pressiodemoapps::Euler2d::Riemann:{
	if( m_icIdentifier == 1){
	  riemann2dIC1(initialState, m_meshObj, m_gamma);
	  return initialState;
	}
	else if (m_icIdentifier == 2){
	  riemann2dIC2(initialState, m_meshObj, m_gamma);
	  return initialState;
	}
	else{
	  throw std::runtime_error("Euler2d Riemann: invalid IC identifier");
	}
      }

      case ::pressiodemoapps::Euler2d::NormalShock:{
	normalShock2dIC(initialState, m_meshObj, m_gamma);
	return initialState;
      }

      case ::pressiodemoapps::Euler2d::DoubleMachReflection:{
	doubleMachReflection2dIC(initialState, m_meshObj, m_gamma);
	return initialState;
      }

      case ::pressiodemoapps::Euler2d::testingonlyneumann:{
	return initialState;
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
	for (int i=1; i<=4; ++i)
	{
	  const auto gID = graph(cell, i);

	  if (gID != -1){
	    const auto ci = gID*numDofPerCell;
	    for (int k=0; k<numDofPerCell; ++k){
	      for (int j=0; j<numDofPerCell; ++j){
		trList.push_back( Tr(jacRowOfCurrentCellDensity+k, ci+j, val0) );
	      }
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
    if (m_probEn == ::pressiodemoapps::Euler2d::SedovFull or
	m_probEn == ::pressiodemoapps::Euler2d::Riemann or
	m_probEn == ::pressiodemoapps::Euler2d::testingonlyneumann)
    {
      using ghost_filler_t  = ::pressiodemoapps::impl::Ghost2dNeumannFiller<
	numDofPerCell, U_t, MeshType, ghost_container_type>;
      ghost_filler_t ghF(stencilSize, U, m_meshObj,
			 m_ghostLeft, m_ghostFront,
			 m_ghostRight, m_ghostBack);

      const auto & rowsBd = m_meshObj.graphRowsOfCellsNearBd();
      for (int it=0; it<rowsBd.size(); ++it){
	ghF(rowsBd[it], it);
      }
    }

    else if (m_probEn == ::pressiodemoapps::Euler2d::SedovSymmetry)
    {
      using ghost_filler_t  = ::pressiodemoapps::impl::Sedov2dSymmetryGhostFiller<
	U_t, MeshType, ghost_container_type>;
      ghost_filler_t ghF(stencilSize, U, m_meshObj,
			 m_ghostLeft, m_ghostFront,
			 m_ghostRight, m_ghostBack);

      const auto & rowsBd = m_meshObj.graphRowsOfCellsNearBd();
      for (int it=0; it<rowsBd.size(); ++it){
	ghF(rowsBd[it], it);
      }
    }

    else if (m_probEn == ::pressiodemoapps::Euler2d::NormalShock)
    {
      using ghost_filler_t =
	::pressiodemoapps::impl::NormalShock2dGhostFiller<
	U_t, MeshType, ghost_container_type>;
      ghost_filler_t ghF(stencilSize, U,
			 currentTime, m_gamma, m_meshObj,
			 m_ghostLeft, m_ghostFront,
			 m_ghostRight, m_ghostBack);

      const auto & rowsBd = m_meshObj.graphRowsOfCellsNearBd();
      for (int it=0; it<rowsBd.size(); ++it){
	ghF(rowsBd[it], it);
      }
    }

    else if (m_probEn == ::pressiodemoapps::Euler2d::DoubleMachReflection)
    {
      using ghost_filler_t =
	::pressiodemoapps::impl::DoubleMachReflection2dGhostFiller<
	U_t, MeshType, ghost_container_type>;
      ghost_filler_t ghF(stencilSize, U,
			 currentTime, m_gamma, m_meshObj,
			 m_ghostLeft, m_ghostFront,
			 m_ghostRight, m_ghostBack);

      const auto & rowsBd = m_meshObj.graphRowsOfCellsNearBd();
      for (int it=0; it<rowsBd.size(); ++it){
	ghF(rowsBd[it], it);
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

    flux_jac_type JLneg, JLpos;
    flux_jac_type JRneg, JRpos;
    flux_jac_type JFneg, JFpos;
    flux_jac_type JBneg, JBpos;

    // flux at faces
    flux_type FLeft(numDofPerCell);
    flux_type FRight(numDofPerCell);
    flux_type FFront(numDofPerCell);
    flux_type FBack(numDofPerCell);

    // reconstructed values at faces
    edge_rec_type uMinusHalfNeg(numDofPerCell);
    edge_rec_type uMinusHalfPos(numDofPerCell);
    edge_rec_type uPlusHalfNeg(numDofPerCell);
    edge_rec_type uPlusHalfPos(numDofPerCell);

    fillGhosts(U, currentTime);

    if (J){
      ::pressiodemoapps::set_zero(*J);
    }

    velocityAndJacInnerCellsImpl(U, currentTime, V, J,
				 FLeft, FRight, FBack, FFront,
				 JLneg, JLpos,
				 JRneg, JRpos,
				 JBneg, JBpos,
				 JFneg, JFpos,
				 uMinusHalfNeg, uMinusHalfPos,
				 uPlusHalfNeg,  uPlusHalfPos);

    velocityAndJacNearBDCellsImpl(U, currentTime, V, J,
				  FLeft, FRight, FBack, FFront,
				  JLneg, JLpos,
				  JRneg, JRpos,
				  JBneg, JBpos,
				  JFneg, JFpos,
				  uMinusHalfNeg, uMinusHalfPos,
				  uPlusHalfNeg,  uPlusHalfPos);
  }

  int findCellBdType(index_t graphRow, int axis) const
  {
    // 0: Neumann
    // 1: Reflective
    // 2: Dirichlet
    constexpr int neumann = 0;
    constexpr int reflective = 1;
    constexpr int dirichlet = 2;

    // we need to change this, since we know which cells are what in advance
    // we can just fit all this logic by having separate loops over the
    // various boundary cells

    if (m_probEn == ::pressiodemoapps::Euler2d::SedovFull or
	m_probEn == ::pressiodemoapps::Euler2d::Riemann or
	m_probEn == ::pressiodemoapps::Euler2d::testingonlyneumann)
    {
      (void)graphRow;
      return neumann;
    }

    else if (m_probEn == ::pressiodemoapps::Euler2d::SedovSymmetry)
    {
      if (axis == 1 && m_meshObj.hasBdLeft2d(graphRow)){
	return reflective;
      }

      if (axis == 2 && m_meshObj.hasBdBack2d(graphRow)){
	return reflective;
      }

      return neumann;
    }

    else if (m_probEn == ::pressiodemoapps::Euler2d::NormalShock)
    {
      if (axis == 2){
	return reflective;
      }
      else{
	return neumann;
      }
    }

    else if (m_probEn == ::pressiodemoapps::Euler2d::DoubleMachReflection)
    {
      const scalar_type one = static_cast<scalar_type>(1);
      const scalar_type six = static_cast<scalar_type>(6);
      const scalar_type wedgePosition = one/six;
      const auto & x = m_meshObj.viewX();
      const auto cellGID = m_meshObj.graph()(graphRow, 0);
      const auto myX = x(cellGID);

      if (axis == 1){
	return neumann;
      }

      if (axis == 2 && m_meshObj.hasBdBack2d(graphRow) && (myX < wedgePosition)){
	return neumann;
      }

      if (axis == 2 && m_meshObj.hasBdBack2d(graphRow) && (myX >= wedgePosition)){
	return reflective;
      }

      return dirichlet;
    }

    return 0;
  }

  template<class U_t, class V_t>
  void velocityAndJacInnerCellsImpl(const U_t & U,
				    const scalar_type currentTime,
				    V_t & V,
				    jacobian_type * J,
				    flux_type & FLeft,
				    flux_type & FRight,
				    flux_type & FBack,
				    flux_type & FFront,
				    flux_jac_type & JLneg,
				    flux_jac_type & JLpos,
				    flux_jac_type & JRneg,
				    flux_jac_type & JRpos,
				    flux_jac_type & JBneg,
				    flux_jac_type & JBpos,
				    flux_jac_type & JFneg,
				    flux_jac_type & JFpos,
				    edge_rec_type & uMinusHalfNeg,
				    edge_rec_type & uMinusHalfPos,
				    edge_rec_type & uPlusHalfNeg,
				    edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;

    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    const auto dxInv   = m_meshObj.dxInv();
    const auto dyInv   = m_meshObj.dyInv();

    using reconstructor_functor_t = pda::impl::ReconstructorFromState<
      dimensionality, edge_rec_type, U_t, MeshType>;

    using jac_fnct_t = pda::impl::FirstOrderInnerCellJacobianFunctor<
      dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, MeshType>;

    // reconstruct functor for face fluxes
    // here we need to use whatever order (m_recEn) user decides
    reconstructor_functor_t ReconstructorX(xAxis, m_recEn, U, m_meshObj,
					   uMinusHalfNeg, uMinusHalfPos,
					   uPlusHalfNeg,  uPlusHalfPos);

    reconstructor_functor_t ReconstructorY(yAxis, m_recEn, U, m_meshObj,
					   uMinusHalfNeg, uMinusHalfPos,
					   uPlusHalfNeg,  uPlusHalfPos);

    // cell jacobian functor
    // currently, REGARDLESS of the reconstruction scheme,
    // we have only first-order Jacobian so we need a first-order reconstructor for the Jacobian
    const auto firstOrder = pda::InviscidFluxReconstruction::FirstOrder;
    edge_rec_type uMinusHalfNegForJ(numDofPerCell);
    edge_rec_type uMinusHalfPosForJ(numDofPerCell);
    edge_rec_type uPlusHalfNegForJ(numDofPerCell);
    edge_rec_type uPlusHalfPosForJ(numDofPerCell);
    reconstructor_functor_t ReconstructorXForJ(xAxis, firstOrder, U, m_meshObj,
					       uMinusHalfNegForJ, uMinusHalfPosForJ,
					       uPlusHalfNegForJ,  uPlusHalfPosForJ);

    reconstructor_functor_t ReconstructorYForJ(yAxis, firstOrder, U, m_meshObj,
					       uMinusHalfNegForJ, uMinusHalfPosForJ,
					       uPlusHalfNegForJ,  uPlusHalfPosForJ);

    jac_fnct_t CellJacobianFunctorX(*J, m_meshObj, JLneg, JLpos, JRneg, JRpos, xAxis);
    jac_fnct_t CellJacobianFunctorY(*J, m_meshObj, JBneg, JBpos, JFneg, JFpos, yAxis);

    const auto & targetGraphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (int it=0; it<targetGraphRows.size(); ++it)
    {
      const auto smPt = targetGraphRows[it];
      const auto vIndex = smPt*numDofPerCell;

      // *** X ***
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
	case pda::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFourDof(FLeft,  uMinusHalfNeg, uMinusHalfPos, normalX_, m_gamma);
	  eeRusanovFluxFourDof(FRight, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gamma);

	  if (J){
	    eeRusanovFluxJacobianFourDof(JLneg, JLpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
					 normalX_, m_gamma);
	    eeRusanovFluxJacobianFourDof(JRneg, JRpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
					 normalX_, m_gamma);
	  }
	  break;
	}

      if (J){
	CellJacobianFunctorX(smPt);
      }

      // *** Y ***
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
	case pda::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFourDof(FBack,  uMinusHalfNeg, uMinusHalfPos, normalY_, m_gamma);
	  eeRusanovFluxFourDof(FFront, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gamma);

	  if (J){
	    eeRusanovFluxJacobianFourDof(JBneg, JBpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
					 normalY_, m_gamma);
	    eeRusanovFluxJacobianFourDof(JFneg, JFpos, uPlusHalfNegForJ,  uPlusHalfPosForJ,
					 normalY_, m_gamma);
	  }
	  break;
	}

      if (J){
	CellJacobianFunctorY(smPt);
      }

      V(vIndex)   = dxInv*(FLeft(0) - FRight(0)) + dyInv*(FBack(0) - FFront(0));
      V(vIndex+1) = dxInv*(FLeft(1) - FRight(1)) + dyInv*(FBack(1) - FFront(1));
      V(vIndex+2) = dxInv*(FLeft(2) - FRight(2)) + dyInv*(FBack(2) - FFront(2));
      V(vIndex+3) = dxInv*(FLeft(3) - FRight(3)) + dyInv*(FBack(3) - FFront(3));
    }
  }

  template<class U_t, class V_t>
  void velocityAndJacNearBDCellsImpl(const U_t & U,
				     const scalar_type currentTime,
				     V_t & V,
				     jacobian_type * J,
				     flux_type & FLeft,
				     flux_type & FRight,
				     flux_type & FBack,
				     flux_type & FFront,
				     flux_jac_type & JLneg,
				     flux_jac_type & JLpos,
				     flux_jac_type & JRneg,
				     flux_jac_type & JRpos,
				     flux_jac_type & JBneg,
				     flux_jac_type & JBpos,
				     flux_jac_type & JFneg,
				     flux_jac_type & JFpos,
				     edge_rec_type & uMinusHalfNeg,
				     edge_rec_type & uMinusHalfPos,
				     edge_rec_type & uPlusHalfNeg,
				     edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;

    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    const auto dxInv   = m_meshObj.dxInv();
    const auto dyInv   = m_meshObj.dyInv();

    using stencil_filler_t = pda::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_container_type,
      U_t, MeshType, ghost_container_type>;
    using rec_fnct_t = pda::impl::ReconstructorFromStencilFourDofPerCell
      <edge_rec_type, stencil_container_type>;

    using jac_fnct_t = pda::impl::FirstOrderBdCellJacobianFunctor<
      dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, MeshType>;

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
    const auto cellJacOrdEn = pda::InviscidFluxReconstruction::FirstOrder;

    stencil_container_type stencilValsForJ = {};
    const auto stencilSizeForJ = reconstructionTypeToStencilSize(cellJacOrdEn);
    pda::resize(stencilValsForJ, numDofPerCell*stencilSizeForJ);
    stencil_filler_t FillStencilValuesXFunctorForJ(stencilSizeForJ, U, m_meshObj,
						   m_ghostLeft, m_ghostRight,
						   stencilValsForJ, xAxis);
    stencil_filler_t FillStencilValuesYFunctorForJ(stencilSizeForJ, U, m_meshObj,
						   m_ghostBack, m_ghostFront,
						   stencilValsForJ, yAxis);

    edge_rec_type uMinusHalfNegForJ(numDofPerCell);
    edge_rec_type uMinusHalfPosForJ(numDofPerCell);
    edge_rec_type uPlusHalfNegForJ(numDofPerCell);
    edge_rec_type uPlusHalfPosForJ(numDofPerCell);
    rec_fnct_t FaceValuesReconstructFunctorForJ(cellJacOrdEn, stencilValsForJ,
						uMinusHalfNegForJ, uMinusHalfPosForJ,
						uPlusHalfNegForJ,  uPlusHalfPosForJ);

    jac_fnct_t CellJacobianFunctorX(*J, m_meshObj, JLneg, JLpos, JRneg, JRpos, xAxis);
    jac_fnct_t CellJacobianFunctorY(*J, m_meshObj, JBneg, JBpos, JFneg, JFpos, yAxis);

    std::array<scalar_type, numDofPerCell> bcCellJacFactorsDefault;
    std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveX;
    std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveY;
    bcCellJacFactorsDefault.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsReflectiveX.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsReflectiveY.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsReflectiveX[1] = static_cast<scalar_type>(-1);
    bcCellJacFactorsReflectiveY[2] = static_cast<scalar_type>(-1);

    // ----------------------------------------
    // loop over cells
    // ----------------------------------------
    const auto & graph = m_meshObj.graph();
    const auto & targetRows = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<targetRows.size(); ++it)
    {
      const auto smPt = targetRows[it];
      const auto vIndex = smPt*numDofPerCell;

      // *** X ***
      StencilFillerX(smPt, it);
      Reconstructor();

      if(J){
	FillStencilValuesXFunctorForJ(smPt, it);
	FaceValuesReconstructFunctorForJ();
      }

      switch(m_fluxEn)
	{
	case pda::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFourDof(FLeft, uMinusHalfNeg, uMinusHalfPos, normalX_, m_gamma);
	  eeRusanovFluxFourDof(FRight, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gamma);

	  if(J){
	    eeRusanovFluxJacobianFourDof(JLneg, JLpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
					 normalX_, m_gamma);
	    eeRusanovFluxJacobianFourDof(JRneg, JRpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
					 normalX_, m_gamma);
	  }
	  break;
	}

      if(J){
	auto bcType = findCellBdType(smPt, xAxis);
	const auto & factorsX = (bcType == 1) ?
	  bcCellJacFactorsReflectiveX : bcCellJacFactorsDefault;
	CellJacobianFunctorX(smPt, factorsX, bcType);
      }

      // *** Y ***
      StencilFillerY(smPt, it);
      Reconstructor();

      if(J){
	FillStencilValuesYFunctorForJ(smPt, it);
	FaceValuesReconstructFunctorForJ();
      }

      switch(m_fluxEn)
	{
	case pda::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFourDof(FBack, uMinusHalfNeg, uMinusHalfPos, normalY_, m_gamma);
	  eeRusanovFluxFourDof(FFront, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gamma);

	  if(J){
	    eeRusanovFluxJacobianFourDof(JBneg, JBpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
					 normalY_, m_gamma);
	    eeRusanovFluxJacobianFourDof(JFneg, JFpos, uPlusHalfNegForJ,  uPlusHalfPosForJ,
					 normalY_, m_gamma);
	  }
	  break;
	}

      if(J){
	auto bcType = findCellBdType(smPt, yAxis);
	const auto & factorsY = (bcType == 1) ?
	  bcCellJacFactorsReflectiveY : bcCellJacFactorsDefault;
	CellJacobianFunctorY(smPt, factorsY, bcType);
      }

      // fill velocity
      V(vIndex)   = dxInv*(FLeft(0) - FRight(0)) + dyInv*(FBack(0) - FFront(0));
      V(vIndex+1) = dxInv*(FLeft(1) - FRight(1)) + dyInv*(FBack(1) - FFront(1));
      V(vIndex+2) = dxInv*(FLeft(2) - FRight(2)) + dyInv*(FBack(2) - FFront(2));
      V(vIndex+3) = dxInv*(FLeft(3) - FRight(3)) + dyInv*(FBack(3) - FFront(3));
    }
  }

  void allocateStencilValuesContainer()
  {
    // the stencil size needed is determined by the desired reconstruction
    // kind NOT from the mesh. THis is important because for example
    // the mesh can have a wider connectivity that what is needed.
    const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_recEn);
    ::pressiodemoapps::resize(m_stencilVals, numDofPerCell*stencilSizeNeeded);
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

    const index_t s1 = m_meshObj.numCellsBd();
    ::pressiodemoapps::resize(m_ghostLeft,  s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostFront, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostBack,  s1, numGhostValues);
  }

protected:
  scalar_type m_gamma = static_cast<scalar_type>(1.4);

  const MeshType & m_meshObj;
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

  mutable stencil_container_type m_stencilVals;
  mutable ghost_container_type m_ghostLeft;
  mutable ghost_container_type m_ghostFront;
  mutable ghost_container_type m_ghostRight;
  mutable ghost_container_type m_ghostBack;

  const std::array<scalar_type, 2> normalX_{1, 0};
  const std::array<scalar_type, 2> normalY_{0, 1};
};

}}}//end namespace
#endif
