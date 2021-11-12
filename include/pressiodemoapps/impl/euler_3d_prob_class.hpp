
#ifndef PRESSIODEMOAPPS_EULER3D_HPP_
#define PRESSIODEMOAPPS_EULER3D_HPP_

// note that the code below is pretty ugly for now, but it works.
// this is inside impl namespace for a reason, and will need to be
// improved later on but we have a starting point.

#include "euler_fluxes.hpp"
#include "euler_flux_jacobian.hpp"
#include "euler_rankine_hugoniot.hpp"
#include "euler_3d_initial_condition.hpp"
#include "functor_ghost_fill_neumann.hpp"
#include "euler_3d_ghost_filler.hpp"
#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "cell_jacobian_first_order.hpp"

namespace pressiodemoapps{ namespace ee{ namespace impl{

template<
  class scalar_t,
  class mesh_t,
  class state_t,
  class velo_t,
  class ghost_t,
  class jacobian_t = void
  >
class Euler3dAppT
{

public:
  using index_t		 = typename mesh_t::index_t;
  using scalar_type	 = scalar_t;
  using state_type	 = state_t;
  using velocity_type	 = velo_t;
  using ghost_type	 = ghost_t;
  using stencil_values_t = state_type;
  using flux_t		 = state_type;
  using edge_rec_t	 = state_type;

  static constexpr int dimensionality{3};
  static constexpr index_t numDofPerCell{5};

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
  using jacobian_type	 = jacobian_t;
  using flux_jac_type = Eigen::Matrix<scalar_type, numDofPerCell, numDofPerCell>;
#endif

public:
#if not defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  Euler3dAppT(const mesh_t & meshObj,
	      ::pressiodemoapps::Euler3d probEn,
	      ::pressiodemoapps::InviscidFluxReconstruction recEn,
	      ::pressiodemoapps::InviscidFluxScheme fluxEnum,
	      int icIdentifier = 1)
    : m_probEn(probEn),
      m_recEn(recEn),
      m_fluxEn(fluxEnum),
      m_icIdentifier(icIdentifier),
      m_meshObj(meshObj)
  {
    computeDofs();
    allocateStencilValuesContainer();
    allocateGhosts();

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
    m_jacobian.resize(m_numDofSampleMesh, m_numDofStencilMesh);
    initializeJacobian();
#endif
  }

#else
  // note that when doing bindings, I need to first construct
  // ghost with {1,1} just so that the numpy array picks up they
  // are 2dim array otherwise it thinks they are 1d array.
  // The right allocation for these is then done inside allocateGhosts.
  Euler3dAppT(const mesh_t & meshObj,
	      ::pressiodemoapps::Euler3d probEn,
	      ::pressiodemoapps::InviscidFluxReconstruction recEn,
	      ::pressiodemoapps::InviscidFluxScheme fluxEnum,
	      int icIdentifier)
    : m_probEn(probEn),
      m_recEn(recEn),
      m_fluxEn(fluxEnum),
      m_icIdentifier(icIdentifier),
      m_meshObj(meshObj),
      m_stencilVals(1),
      m_ghostLeft({1,1}),
      m_ghostRight({1,1}),
      m_ghostBack({1,1}),
      m_ghostFront({1,1}),
      m_ghostBottom({1,1}),
      m_ghostTop({1,1})
  {
    computeDofs();
    allocateStencilValuesContainer();
    allocateGhosts();
  }
#endif

  state_type initialCondition() const
  {
    state_type IC(m_numDofStencilMesh);

    switch(m_probEn)
      {
      case ::pressiodemoapps::Euler3d::PeriodicSmooth:{
	euler3dsmoothInitialCondition(IC, m_meshObj, m_gamma);
	return IC;
      }
      case ::pressiodemoapps::Euler3d::SedovSymmetry:{
	sedov3dInitialCondition(IC, m_meshObj, m_gamma);
	return IC;
      }
      };

    return IC;
  }

  scalar_type gamma()           const{ return m_gamma; }
  index_t totalDofSampleMesh()  const{ return m_numDofSampleMesh; }
  index_t totalDofStencilMesh() const{ return m_numDofStencilMesh; }

  velocity_type createVelocity() const
  {
    velocity_type V(m_numDofSampleMesh);
    return V;
  }

  void velocity(const state_type & state,
		const scalar_type timeValue,
		velocity_type & V) const
  {
    flux_t FL(numDofPerCell);
    flux_t FR(numDofPerCell);
    flux_t FB(numDofPerCell);
    flux_t FF(numDofPerCell);
    flux_t FD(numDofPerCell);
    flux_t FU(numDofPerCell);
    edge_rec_t uMinusHalfNeg(numDofPerCell);
    edge_rec_t uMinusHalfPos(numDofPerCell);
    edge_rec_t uPlusHalfNeg (numDofPerCell);
    edge_rec_t uPlusHalfPos (numDofPerCell);

    fillGhosts(state);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
    if (!m_onlyComputeVelocity)
    {
      auto values = m_jacobian.valuePtr();
      for (int i=0; i<m_jacobian.nonZeros(); ++i){
	values[i] = 0;
      }
    }
#endif

    velocityCellsNearBdImpl(state, timeValue,
			    V, FL, FR, FB, FF, FD, FU,
			    uMinusHalfNeg, uMinusHalfPos,
			    uPlusHalfNeg,  uPlusHalfPos);
    velocityInnerCellsImpl(state, timeValue,
			   V, FL, FR, FB, FF, FD, FU,
			   uMinusHalfNeg, uMinusHalfPos,
			   uPlusHalfNeg,  uPlusHalfPos);
  }

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
  jacobian_type createJacobian() const{
    jacobian_type JJ(m_numDofSampleMesh, m_numDofStencilMesh);
    return JJ;
  }

  void jacobian(const state_type & state,
		const scalar_type time,
		jacobian_type & J) const
  {
    if (!m_onlyComputeVelocity){
      // relies on jacobian been computed in velocity
      J = m_jacobian;
    }
  }

  // the Jacobian is by default fused with the velocity,
  // this method allows one to disable the jacobian
  // so only velocity is computed
  void disableJacobian() {
    m_onlyComputeVelocity = true;
  }
#endif

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

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
  void initializeJacobian()
  {
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
	const auto D0 = graph(cell, 5);
	const auto U0 = graph(cell, 6);

	for (int i=1; i<=6; ++i)
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

    m_jacobian.setFromTriplets(trList.begin(), trList.end());
  }
#endif

  void velocityInnerCellsImpl(const state_type & U,
			      const scalar_type timeIn,
			      velocity_type & V,
			      flux_t & FL,
			      flux_t & FR,
			      flux_t & FB,
			      flux_t & FF,
			      flux_t & FD,
			      flux_t & FU,
			      edge_rec_t & uMinusHalfNeg,
			      edge_rec_t & uMinusHalfPos,
			      edge_rec_t & uPlusHalfNeg,
			      edge_rec_t & uPlusHalfPos) const
  {

    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    constexpr int zAxis = 3;
    const auto dxInv   = m_meshObj.dxInv();
    const auto dyInv   = m_meshObj.dyInv();
    const auto dzInv   = m_meshObj.dzInv();

    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromState<
      dimensionality, edge_rec_t, state_type, mesh_t>;

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
    using jac_fnct_t = ::pressiodemoapps::impl::FirstOrderInnerCellJacobianFunctor<
      dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, mesh_t>;
#endif

    // reconstruct functor for face fluxes
    // here we need to use whatever order (m_recEn) user decides
    rec_fnct_t ReconstructorX(xAxis, m_recEn, U, m_meshObj,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);

    rec_fnct_t ReconstructorY(yAxis, m_recEn, U, m_meshObj,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);

    rec_fnct_t ReconstructorZ(zAxis, m_recEn, U, m_meshObj,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
    // cell jacobian functor
    // currently, REGARDLESS of the reconstruction scheme,
    // we have only first-order Jacobian so we need a first-order reconstructor for the Jacobian
    flux_jac_type JLneg, JLpos;
    flux_jac_type JRneg, JRpos;
    flux_jac_type JBneg, JBpos;
    flux_jac_type JFneg, JFpos;
    flux_jac_type JDneg, JDpos;
    flux_jac_type JUneg, JUpos;

    const auto firstOrder = ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder;
    edge_rec_t uMinusHalfNegForJ(numDofPerCell);
    edge_rec_t uMinusHalfPosForJ(numDofPerCell);
    edge_rec_t uPlusHalfNegForJ(numDofPerCell);
    edge_rec_t uPlusHalfPosForJ(numDofPerCell);
    rec_fnct_t ReconstructorXForJ(xAxis, firstOrder, U, m_meshObj,
				  uMinusHalfNegForJ, uMinusHalfPosForJ,
				  uPlusHalfNegForJ,  uPlusHalfPosForJ);

    rec_fnct_t ReconstructorYForJ(yAxis, firstOrder, U, m_meshObj,
				  uMinusHalfNegForJ, uMinusHalfPosForJ,
				  uPlusHalfNegForJ,  uPlusHalfPosForJ);

    rec_fnct_t ReconstructorZForJ(zAxis, firstOrder, U, m_meshObj,
				  uMinusHalfNegForJ, uMinusHalfPosForJ,
				  uPlusHalfNegForJ,  uPlusHalfPosForJ);

    jac_fnct_t CellJacobianFunctorX(m_jacobian, m_meshObj, JLneg, JLpos, JRneg, JRpos, xAxis);
    jac_fnct_t CellJacobianFunctorY(m_jacobian, m_meshObj, JBneg, JBpos, JFneg, JFpos, yAxis);
    jac_fnct_t CellJacobianFunctorZ(m_jacobian, m_meshObj, JDneg, JDpos, JUneg, JUpos, zAxis);
#endif

    // deal with cells away from boundaries
    const auto & rowsIn = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (int it=0; it<rowsIn.size(); ++it)
    {
      const auto smPt = rowsIn[it];
      const auto vIndex = smPt*numDofPerCell;

      // *** X ***
      ReconstructorX.template operator()<numDofPerCell>(smPt);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (!m_onlyComputeVelocity){
	// note that REGARDLESS of the reconstruction scheme,
	// we currently only have only first-order Jacobian so we need
	// to run the reconstructor for the Jacobian
	// which will ensure that uMinusNegForJ, etc have the right values
	ReconstructorXForJ.template operator()<numDofPerCell>(smPt);
      }
#endif

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFiveDof(FL, uMinusHalfNeg, uMinusHalfPos, normalX_, m_gamma);
	  eeRusanovFluxFiveDof(FR, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gamma);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
	  if (!m_onlyComputeVelocity){
	    eeRusanovFluxJacobianFiveDof(JLneg, JLpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
					 normalX_, m_gamma);
	    eeRusanovFluxJacobianFiveDof(JRneg, JRpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
					 normalX_, m_gamma);
	  }
#endif
	  break;
	}

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (!m_onlyComputeVelocity){
	CellJacobianFunctorX(smPt);
      }
#endif

      // *** Y ***
      ReconstructorY.template operator()<numDofPerCell>(smPt);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (!m_onlyComputeVelocity){
	// note that REGARDLESS of the reconstruction scheme,
	// we currently only have only first-order Jacobian so we need
	// to run the reconstructor for the Jacobian
	// which will ensure that uMinusNegForJ, etc have the right values
	ReconstructorYForJ.template operator()<numDofPerCell>(smPt);
      }
#endif

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFiveDof(FB, uMinusHalfNeg, uMinusHalfPos, normalY_, m_gamma);
	  eeRusanovFluxFiveDof(FF, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gamma);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
	  if (!m_onlyComputeVelocity){
	    eeRusanovFluxJacobianFiveDof(JBneg, JBpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
					 normalY_, m_gamma);
	    eeRusanovFluxJacobianFiveDof(JFneg, JFpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
					 normalY_, m_gamma);
	  }
#endif
	  break;
	}

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (!m_onlyComputeVelocity){
	CellJacobianFunctorY(smPt);
      }
#endif

      // *** Z ***
      ReconstructorZ.template operator()<numDofPerCell>(smPt);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (!m_onlyComputeVelocity){
	// note that REGARDLESS of the reconstruction scheme,
	// we currently only have only first-order Jacobian so we need
	// to run the reconstructor for the Jacobian
	// which will ensure that uMinusNegForJ, etc have the right values
	ReconstructorZForJ.template operator()<numDofPerCell>(smPt);
      }
#endif

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFiveDof(FD, uMinusHalfNeg, uMinusHalfPos, normalZ_, m_gamma);
	  eeRusanovFluxFiveDof(FU, uPlusHalfNeg,  uPlusHalfPos,  normalZ_, m_gamma);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
	  if (!m_onlyComputeVelocity){
	    eeRusanovFluxJacobianFiveDof(JDneg, JDpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
					 normalZ_, m_gamma);
	    eeRusanovFluxJacobianFiveDof(JUneg, JUpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
					 normalZ_, m_gamma);
	  }
#endif
	  break;
	}

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (!m_onlyComputeVelocity){
	CellJacobianFunctorZ(smPt);
      }
#endif

      V(vIndex)   = dxInv*(FL(0) - FR(0)) + dyInv*(FB(0) - FF(0)) + dzInv*(FD(0) - FU(0));
      V(vIndex+1) = dxInv*(FL(1) - FR(1)) + dyInv*(FB(1) - FF(1)) + dzInv*(FD(1) - FU(1));
      V(vIndex+2) = dxInv*(FL(2) - FR(2)) + dyInv*(FB(2) - FF(2)) + dzInv*(FD(2) - FU(2));
      V(vIndex+3) = dxInv*(FL(3) - FR(3)) + dyInv*(FB(3) - FF(3)) + dzInv*(FD(3) - FU(3));
      V(vIndex+4) = dxInv*(FL(4) - FR(4)) + dyInv*(FB(4) - FF(4)) + dzInv*(FD(4) - FU(4));
    }
  }

  int findCellBdType(index_t graphRow, int axis) const
  {
    // 0: Neumann
    // 1: Reflective
    // 2: Dirichlet
    constexpr int neumann = 0;
    constexpr int reflective = 1;
    constexpr int dirichlet = 2;

    if (m_probEn == ::pressiodemoapps::Euler3d::SedovSymmetry)
    {
      if (axis == 1 && m_meshObj.hasBdLeft3d(graphRow)){
	return reflective;
      }

      if (axis == 2 && m_meshObj.hasBdBack3d(graphRow)){
	return reflective;
      }

      if (axis == 3 && m_meshObj.hasBdBottom3d(graphRow)){
	return reflective;
      }

      return neumann;
    }

    return 0;
  }

  void velocityCellsNearBdImpl(const state_type & U,
			       const scalar_type timeIn,
			       velocity_type & V,
			       flux_t & FL,
			       flux_t & FR,
			       flux_t & FB,
			       flux_t & FF,
			       flux_t & FD,
			       flux_t & FU,
			       edge_rec_t & uMinusHalfNeg,
			       edge_rec_t & uMinusHalfPos,
			       edge_rec_t & uPlusHalfNeg,
			       edge_rec_t & uPlusHalfPos) const
  {

    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    constexpr int zAxis = 3;

    const auto dxInv   = m_meshObj.dxInv();
    const auto dyInv   = m_meshObj.dyInv();
    const auto dzInv   = m_meshObj.dzInv();

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
    flux_jac_type JLneg, JLpos;
    flux_jac_type JRneg, JRpos;
    flux_jac_type JBneg, JBpos;
    flux_jac_type JFneg, JFpos;
    flux_jac_type JDneg, JDpos;
    flux_jac_type JUneg, JUpos;
#endif

    using sfiller_t = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_values_t, state_type, mesh_t, ghost_t>;

    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromStencilFiveDofPerCell<
      edge_rec_t, stencil_values_t>;

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
    using jac_fnct_t = ::pressiodemoapps::impl::FirstOrderBdCellJacobianFunctor<
      dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, mesh_t>;
#endif

    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    sfiller_t StencilFillerX(stencilSize, U, m_meshObj,
			     m_ghostLeft, m_ghostRight,
			     m_stencilVals, xAxis);

    sfiller_t StencilFillerY(stencilSize, U, m_meshObj,
			     m_ghostBack, m_ghostFront,
			     m_stencilVals, yAxis);

    sfiller_t StencilFillerZ(stencilSize, U, m_meshObj,
			     m_ghostBottom, m_ghostTop,
			     m_stencilVals, zAxis);

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
    sfiller_t FillStencilValuesXFunctorForJ(stencilSizeForJ, U, m_meshObj,
					    m_ghostLeft, m_ghostRight, stencilValsForJ, xAxis);

    sfiller_t FillStencilValuesYFunctorForJ(stencilSizeForJ, U, m_meshObj,
					    m_ghostBack, m_ghostFront, stencilValsForJ, yAxis);

    sfiller_t FillStencilValuesZFunctorForJ(stencilSizeForJ, U, m_meshObj,
					    m_ghostBottom, m_ghostTop, stencilValsForJ, zAxis);

    edge_rec_t uMinusHalfNegForJ(numDofPerCell);
    edge_rec_t uMinusHalfPosForJ(numDofPerCell);
    edge_rec_t uPlusHalfNegForJ(numDofPerCell);
    edge_rec_t uPlusHalfPosForJ(numDofPerCell);
    rec_fnct_t FaceValuesReconstructFunctorForJ(cellJacOrdEn, stencilValsForJ,
						uMinusHalfNegForJ, uMinusHalfPosForJ,
						uPlusHalfNegForJ,  uPlusHalfPosForJ);

    jac_fnct_t CellJacobianFunctorX(m_jacobian, m_meshObj, JLneg, JLpos, JRneg, JRpos, xAxis);
    jac_fnct_t CellJacobianFunctorY(m_jacobian, m_meshObj, JBneg, JBpos, JFneg, JFpos, yAxis);
    jac_fnct_t CellJacobianFunctorZ(m_jacobian, m_meshObj, JDneg, JDpos, JUneg, JUpos, zAxis);

    std::array<scalar_type, numDofPerCell> bcCellJacFactorsDefault;
    bcCellJacFactorsDefault.fill(static_cast<scalar_type>(1));

    std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveX;
    std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveY;
    std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveZ;

    bcCellJacFactorsReflectiveX.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsReflectiveY.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsReflectiveZ.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsReflectiveX[1] = static_cast<scalar_type>(-1);
    bcCellJacFactorsReflectiveY[2] = static_cast<scalar_type>(-1);
    bcCellJacFactorsReflectiveZ[3] = static_cast<scalar_type>(-1);
#endif

    // -----
    // loop
    // -----
    const auto & specialRows = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<specialRows.size(); ++it)
    {
      const auto smPt = specialRows[it];

      // ------------
      // X
      // ------------
      StencilFillerX(smPt, it);
      Reconstructor();

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (!m_onlyComputeVelocity){
	FillStencilValuesXFunctorForJ(smPt, it);
	FaceValuesReconstructFunctorForJ();
      }
#endif

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFiveDof(FL, uMinusHalfNeg, uMinusHalfPos, normalX_, m_gamma);
	  eeRusanovFluxFiveDof(FR, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gamma);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
	  if (!m_onlyComputeVelocity){
	    eeRusanovFluxJacobianFiveDof(JLneg, JLpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
					 normalX_, m_gamma);
	    eeRusanovFluxJacobianFiveDof(JRneg, JRpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
					 normalX_, m_gamma);
	  }
#endif

	  break;
	}

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (!m_onlyComputeVelocity){
	auto bcType = findCellBdType(smPt, xAxis);
	const auto & factorsX = (bcType == 1) ? bcCellJacFactorsReflectiveX : bcCellJacFactorsDefault;
	CellJacobianFunctorX(smPt, factorsX, bcType);
      }
#endif

      // ------------
      // Y
      // ------------
      StencilFillerY(smPt, it);
      Reconstructor();

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (!m_onlyComputeVelocity){
      FillStencilValuesYFunctorForJ(smPt, it);
      FaceValuesReconstructFunctorForJ();
      }
#endif

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFiveDof(FB, uMinusHalfNeg, uMinusHalfPos, normalY_, m_gamma);
	  eeRusanovFluxFiveDof(FF, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gamma);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
	  if (!m_onlyComputeVelocity){
	  eeRusanovFluxJacobianFiveDof(JBneg, JBpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
				       normalY_, m_gamma);
	  eeRusanovFluxJacobianFiveDof(JFneg, JFpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
				       normalY_, m_gamma);
	  }
#endif

	  break;
	}

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (!m_onlyComputeVelocity){
	auto bcType = findCellBdType(smPt, yAxis);
	const auto & factorsY = (bcType == 1) ? bcCellJacFactorsReflectiveY : bcCellJacFactorsDefault;
	CellJacobianFunctorY(smPt, factorsY, bcType);
      }
#endif

      // ------------
      // Z
      // ------------
      StencilFillerZ(smPt, it);
      Reconstructor();

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (!m_onlyComputeVelocity){
	FillStencilValuesZFunctorForJ(smPt, it);
	FaceValuesReconstructFunctorForJ();
      }
#endif

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFiveDof(FD, uMinusHalfNeg, uMinusHalfPos, normalZ_, m_gamma);
	  eeRusanovFluxFiveDof(FU, uPlusHalfNeg,  uPlusHalfPos,  normalZ_, m_gamma);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
	  if (!m_onlyComputeVelocity){
	  eeRusanovFluxJacobianFiveDof(JDneg, JDpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
				       normalZ_, m_gamma);
	  eeRusanovFluxJacobianFiveDof(JUneg, JUpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
				       normalZ_, m_gamma);
	  }
#endif

	  break;
	}

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (!m_onlyComputeVelocity){
	auto bcType = findCellBdType(smPt, zAxis);
	const auto & factorsZ = (bcType == 1) ? bcCellJacFactorsReflectiveZ : bcCellJacFactorsDefault;
	CellJacobianFunctorZ(smPt, factorsZ, bcType);
      }
#endif

      const auto vIndex = smPt*numDofPerCell;
      V(vIndex)   = dxInv*(FL(0) - FR(0)) + dyInv*(FB(0) - FF(0)) + dzInv*(FD(0) - FU(0));
      V(vIndex+1) = dxInv*(FL(1) - FR(1)) + dyInv*(FB(1) - FF(1)) + dzInv*(FD(1) - FU(1));
      V(vIndex+2) = dxInv*(FL(2) - FR(2)) + dyInv*(FB(2) - FF(2)) + dzInv*(FD(2) - FU(2));
      V(vIndex+3) = dxInv*(FL(3) - FR(3)) + dyInv*(FB(3) - FF(3)) + dzInv*(FD(3) - FU(3));
      V(vIndex+4) = dxInv*(FL(4) - FR(4)) + dyInv*(FB(4) - FF(4)) + dzInv*(FD(4) - FU(4));
    }
  }


  void fillGhosts(const state_type & U) const
  {
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    if (m_probEn == ::pressiodemoapps::Euler3d::SedovSymmetry)
    {
      using ghost_filler_t  = ::pressiodemoapps::ee::impl::Ghost3dSedov<state_type, mesh_t, ghost_t>;
      ghost_filler_t ghF(stencilSize,   U, m_meshObj,
			 m_ghostLeft,   m_ghostRight,
			 m_ghostBack,   m_ghostFront,
			 m_ghostBottom, m_ghostTop);

      const auto & rowsBd = m_meshObj.graphRowsOfCellsNearBd();

      if (stencilSize==3){
	for (int it=0; it<rowsBd.size(); ++it){
	  ghF.template operator()<3>(rowsBd[it], it);
	}
      }
      else if (stencilSize==5){
	for (int it=0; it<rowsBd.size(); ++it){
	  ghF.template operator()<5>(rowsBd[it], it);
	}
      }else{
	throw std::runtime_error("missing impl");
      }
    }
  }

private:
  void allocateGhosts()
  {
    /*
      for stencil = 3, at leftBoundary:
      ---------------
      |	 0,1,2,3,4 ||
      | rho,       ||
      | rho u,	   ||
      | rho v,     ||
      | rho w,     ||
      | E	   ||
      ---------------

      for stencil = 7, at leftBoundary:
      --------------------------------------
      |	10,11,12,13,14  | 5,6,7,8,9 |  0,1,2,3,4 ||
      |	     	        |
      | rho,            | rho,       | rho        ||
      | rho u,	        | rho*u      | rho*u      ||
      | rho v,          | rho*v      | rho*v      ||
      | rho w,          | rho*w      | rho*w      ||
      | E	        | E          | E          ||
      ---------------------------------------
     */

    const auto stencilSize    = reconstructionTypeToStencilSize(m_recEn);
    const auto numGhostValues = numDofPerCell*((stencilSize-1)/2);

    const index_t s1 = m_meshObj.numCellsBd();
    ::pressiodemoapps::resize(m_ghostLeft,   s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight,  s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostBack,   s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostFront,  s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostBottom, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostTop,    s1, numGhostValues);
  }

private:
  const scalar_type m_gamma{1.4};

  ::pressiodemoapps::Euler3d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  ::pressiodemoapps::InviscidFluxScheme m_fluxEn;

  int m_icIdentifier = 1;

  const mesh_t & m_meshObj;
  mutable stencil_values_t m_stencilVals;

  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  mutable ghost_t m_ghostLeft;
  mutable ghost_t m_ghostRight;
  mutable ghost_t m_ghostBack;
  mutable ghost_t m_ghostFront;
  mutable ghost_t m_ghostBottom;
  mutable ghost_t m_ghostTop;

  const std::array<scalar_type, dimensionality> normalX_{1, 0, 0};
  const std::array<scalar_type, dimensionality> normalY_{0, 1, 0};
  const std::array<scalar_type, dimensionality> normalZ_{0, 0, 1};

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
  mutable jacobian_type m_jacobian = {};
  std::size_t m_jacNonZerosCount = {};
  bool m_onlyComputeVelocity = false;
#endif
};

}}}//end namespace
#endif
