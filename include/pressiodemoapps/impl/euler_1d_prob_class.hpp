
#ifndef PRESSIODEMOAPPS_EULER1D_APP_HPP_
#define PRESSIODEMOAPPS_EULER1D_APP_HPP_

#include "euler_fluxes.hpp"
#include "euler_flux_jacobian.hpp"
#include "euler_rankine_hugoniot.hpp"
#include "euler_1d_initial_condition.hpp"
#include "euler_1d_ghost_filler.hpp"
#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "cell_jacobian_first_order.hpp"
#include "Eigen/Sparse"

namespace pressiodemoapps{ namespace implee1d{

template<class MeshType>
class EigenEuler1dApp
{

public:
  using index_t	      = typename MeshType::index_t;
  using scalar_type   = typename MeshType::scalar_t;
  using state_type    = Eigen::Matrix<scalar_type,-1,1>;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, index_t>;

  static constexpr int	   dimensionality{1};
  static constexpr index_t numDofPerCell{3};

private:
  using ghost_container_type   = Eigen::Matrix<scalar_type,-1,-1, Eigen::RowMajor>;
  using stencil_container_type = Eigen::Matrix<scalar_type,-1,1>;
  using flux_type	       = Eigen::Matrix<scalar_type,-1,1>;
  using edge_rec_type	       = Eigen::Matrix<scalar_type,-1,1>;
  using flux_jac_type = Eigen::Matrix<scalar_type, numDofPerCell, numDofPerCell>;

public:
  EigenEuler1dApp(const MeshType & meshObj,
		  ::pressiodemoapps::Euler1d probEnum,
		  ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		  ::pressiodemoapps::InviscidFluxScheme fluxEnum)
    : m_meshObj(meshObj), m_probEn(probEnum), m_recEn(recEnum), m_fluxEn(fluxEnum)
  {

    // calculate total num of dofs on sample and stencil mesh
    m_numDofStencilMesh = m_meshObj.stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * numDofPerCell;

    allocateStencilValuesContainer();
    allocateGhostValues();
  }

  scalar_type gamma() const{
    return m_gamma;
  }

  state_type initialCondition() const
  {
    state_type res(m_numDofStencilMesh);

    if (m_probEn == ::pressiodemoapps::Euler1d::PeriodicSmooth){
      ::pressiodemoapps::impl::euler1dsineInitialCondition(res, m_meshObj, m_gamma);
    }
    else if (m_probEn == ::pressiodemoapps::Euler1d::Sod){
      ::pressiodemoapps::impl::sod1dInitialCondition(res, m_meshObj, m_gamma);
    }
    else if (m_probEn == ::pressiodemoapps::Euler1d::Lax){
      ::pressiodemoapps::impl::lax1dInitialCondition(res, m_meshObj, m_gamma);
    }
    else{
      //nothing
    }
    return res;
  }

protected:
  void initializeJacobian(jacobian_type & J)
  {
    J.resize(m_numDofSampleMesh, m_numDofStencilMesh);
    initializeJacobianFirstOrder(J);

    // compress to make it a real Crs matrix
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
	const auto R0 = graph(cell, 2);
	if (L0 != -1){
	  const auto ciL0 = L0*numDofPerCell;
	  for (int k=0; k<numDofPerCell; ++k){
	    for (int j=0; j<numDofPerCell; ++j){
	      trList.push_back( Tr(jacRowOfCurrentCellDensity+k, ciL0+j, val0) );
	    }
	  }
	}

	if (R0 != -1){
	  const auto ciR0 = R0*numDofPerCell;
	  for (int k=0; k<numDofPerCell; ++k){
	    for (int j=0; j<numDofPerCell; ++j){
	      trList.push_back( Tr(jacRowOfCurrentCellDensity+k, ciR0+j, val0) );
	    }
	  }
	}
      }

    J.setFromTriplets(trList.begin(), trList.end());
  }

  template<class U_t>
  void fillGhostsIfNeeded(const U_t & U,
			  scalar_type /*currTime*/) const
  {
    const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_recEn);
    // only need ghosts for specific problems
    if (m_probEn == ::pressiodemoapps::Euler1d::Sod or
	m_probEn == ::pressiodemoapps::Euler1d::Lax)
    {
      using ghost_filler_t  = ::pressiodemoapps::impl::Ghost1dNeumannFiller<
	numDofPerCell, U_t, MeshType, ghost_container_type>;
      ghost_filler_t ghF(stencilSizeNeeded, U, m_meshObj, m_ghostLeft, m_ghostRight);
      const auto & rowsBd = m_meshObj.graphRowsOfCellsNearBd();
      for (int it=0; it<rowsBd.size(); ++it){
	ghF(rowsBd[it], it);
      }
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

    // reconstructed values at faces
    edge_rec_type uMinusHalfNeg(numDofPerCell);
    edge_rec_type uMinusHalfPos(numDofPerCell);
    edge_rec_type uPlusHalfNeg(numDofPerCell);
    edge_rec_type uPlusHalfPos(numDofPerCell);

    // flux at left and right faces
    flux_type FL(numDofPerCell), FR(numDofPerCell);

    // jac at left face wrt to left,right reconstructed states
    flux_jac_type JLneg, JLpos;
    // jac at right face wrt to left,right reconstructed states
    flux_jac_type JRneg, JRpos;

    fillGhostsIfNeeded(U, currentTime);

    if (J){
      ::pressiodemoapps::set_zero(*J);
    }

    velocityAndJacInnerCellsImpl(U, currentTime,
				 V, J, FL, FR, JLneg, JLpos, JRneg, JRpos,
				 uMinusHalfNeg, uMinusHalfPos,
				 uPlusHalfNeg,  uPlusHalfPos);

    velocityAndJacNearBDCellsImpl(U, currentTime,
				  V, J, FL, FR, JLneg, JLpos, JRneg, JRpos,
				  uMinusHalfNeg, uMinusHalfPos,
				  uPlusHalfNeg,  uPlusHalfPos);
  }

  template<class U_t, class V_t>
  void velocityAndJacInnerCellsImpl(const U_t & U,
				    const scalar_type currentTime,
				    V_t & V,
				    jacobian_type * J,
				    flux_type & FL,
				    flux_type & FR,
				    flux_jac_type & JLneg,
				    flux_jac_type & JLpos,
				    flux_jac_type & JRneg,
				    flux_jac_type & JRpos,
				    edge_rec_type & uMinusHalfNeg,
				    edge_rec_type & uMinusHalfPos,
				    edge_rec_type & uPlusHalfNeg,
				    edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;

    // reconstruct functor for face fluxes
    // here we need to use whatever order (m_recEn) user decides
    using rec_fnct_t = pda::impl::ReconstructorFromState<
      dimensionality, edge_rec_type, U_t, MeshType>;
    rec_fnct_t Reconstructor(m_recEn, U, m_meshObj,
			     uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg,  uPlusHalfPos);

    // cell jacobian functor
    // currently, REGARDLESS of the reconstruction scheme,
    // we have only first-order Jacobian so we need a first-order reconstructor
    // for the Jacobian
    edge_rec_type uMinusHalfNegForJ(numDofPerCell);
    edge_rec_type uMinusHalfPosForJ(numDofPerCell);
    edge_rec_type uPlusHalfNegForJ(numDofPerCell);
    edge_rec_type uPlusHalfPosForJ(numDofPerCell);
    rec_fnct_t ReconstructorForJ(pda::InviscidFluxReconstruction::FirstOrder,
				 U, m_meshObj,
				 uMinusHalfNegForJ, uMinusHalfPosForJ,
				 uPlusHalfNegForJ,  uPlusHalfPosForJ);

    using jac_fnct_t = pda::impl::FirstOrderInnerCellJacobianFunctor<
      dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, MeshType>;
    jac_fnct_t CellJacobianFunctor(*J, m_meshObj, JLneg, JLpos, JRneg, JRpos);

    // loop
    const auto dxInv = m_meshObj.dxInv();
    const auto & targetGraphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (int it=0; it<targetGraphRows.size(); ++it)
      {
	const auto smPt = targetGraphRows[it];
	Reconstructor.template operator()<numDofPerCell>(smPt);

	if(J){
	  // note that REGARDLESS of the reconstruction scheme,
	  // we currently only have only first-order Jacobian so we need
	  // to run the reconstructor for the Jacobian
	  // which will ensure that uMinusNegForJ, etc have the right values
	  ReconstructorForJ.template operator()<numDofPerCell>(smPt);
	}

	switch(m_fluxEn)
	  {
	  case pda::InviscidFluxScheme::Rusanov:
	    ee::impl::eeRusanovFluxThreeDof(FL, uMinusHalfNeg, uMinusHalfPos, m_gamma);
	    ee::impl::eeRusanovFluxThreeDof(FR, uPlusHalfNeg,  uPlusHalfPos,  m_gamma);

	    if(J){
	      ee::impl::eeRusanovFluxJacobianThreeDof(JLneg, JLpos, uMinusHalfNegForJ,
						      uMinusHalfPosForJ, m_gamma);
	      ee::impl::eeRusanovFluxJacobianThreeDof(JRneg, JRpos, uPlusHalfNegForJ,
						      uPlusHalfPosForJ, m_gamma);
	    }
	    break;
	  }

	const auto vIndex = smPt*numDofPerCell;
	V(vIndex)   = dxInv*(FL(0) - FR(0));
	V(vIndex+1) = dxInv*(FL(1) - FR(1));
	V(vIndex+2) = dxInv*(FL(2) - FR(2));

	if(J){
	  CellJacobianFunctor(smPt);
	}
      }
  }

  template<class U_t, class V_t>
  void velocityAndJacNearBDCellsImpl(const U_t & U,
				     const scalar_type currentTime,
				     V_t & V,
				     jacobian_type * J,
				     flux_type & FL,
				     flux_type & FR,
				     flux_jac_type & JLneg,
				     flux_jac_type & JLpos,
				     flux_jac_type & JRneg,
				     flux_jac_type & JRpos,
				     edge_rec_type & uMinusHalfNeg,
				     edge_rec_type & uMinusHalfPos,
				     edge_rec_type & uPlusHalfNeg,
				     edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_container_type,
      U_t, MeshType, ghost_container_type>;

    using rec_functor_t = pda::impl::ReconstructorFromStencilThreeDofPerCell<
      edge_rec_type, stencil_container_type>;

    using jac_functor_t = pda::impl::FirstOrderBdCellJacobianFunctor<
      dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, MeshType>;

    // -----------------------------------------------
    // functors needed to compute face fluxes
    // NOTE: for face fluxes, we must be consistent with whatever scheme
    // was specified by user, so we need to use this enum: m_recEn
    // -----------------------------------------------
    stencil_filler_t FillStencilValuesFunctor(reconstructionTypeToStencilSize(m_recEn), U,
					      m_meshObj, m_ghostLeft, m_ghostRight,
					      m_stencilVals);
    rec_functor_t FaceValuesReconstructFunctor(m_recEn, m_stencilVals,
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
    stencil_filler_t FillStencilValuesFunctorForJ(stencilSizeForJ, U, m_meshObj,
						  m_ghostLeft, m_ghostRight, stencilValsForJ);

    edge_rec_type uMinusHalfNegForJ(numDofPerCell);
    edge_rec_type uMinusHalfPosForJ(numDofPerCell);
    edge_rec_type uPlusHalfNegForJ(numDofPerCell);
    edge_rec_type uPlusHalfPosForJ(numDofPerCell);
    rec_functor_t FaceValuesReconstructFunctorForJ(cellJacOrdEn, stencilValsForJ,
						   uMinusHalfNegForJ, uMinusHalfPosForJ,
						   uPlusHalfNegForJ,  uPlusHalfPosForJ);
    jac_functor_t CellJacobianFunctor(*J, m_meshObj, JLneg, JLpos, JRneg, JRpos);

    std::array<scalar_type, numDofPerCell> bcCellJacFactors;
    bcCellJacFactors.fill(static_cast<scalar_type>(1));

    // ----------------------------------------
    // loop over cells
    // ----------------------------------------
    const auto dxInv = m_meshObj.dxInv();
    const auto & targetGraphRows = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<targetGraphRows.size(); ++it)
    {
      const auto smPt = targetGraphRows[it];

      FillStencilValuesFunctor(smPt, it);
      FaceValuesReconstructFunctor();

      if(J){
	// note that REGARDLESS of the reconstruction scheme,
	// we currently only have only first-order Jacobian so we need
	// to run the reconstructor for the Jacobian
	// which will ensure that uMinusNegForJ, etc have the right values
	FillStencilValuesFunctorForJ(smPt, it);
	FaceValuesReconstructFunctorForJ();
      }

      switch(m_fluxEn)
      {
      case pda::InviscidFluxScheme::Rusanov:
	ee::impl::eeRusanovFluxThreeDof(FL, uMinusHalfNeg, uMinusHalfPos, m_gamma);
	ee::impl::eeRusanovFluxThreeDof(FR, uPlusHalfNeg,  uPlusHalfPos,  m_gamma);

	if (J){
	  ee::impl::eeRusanovFluxJacobianThreeDof(JLneg, JLpos, uMinusHalfNegForJ,
						  uMinusHalfPosForJ, m_gamma);
	  ee::impl::eeRusanovFluxJacobianThreeDof(JRneg, JRpos, uPlusHalfNegForJ,
						  uPlusHalfPosForJ, m_gamma);
	}
	break;
      }

      const auto vIndexCurrentCellDensity = smPt*numDofPerCell;
      V(vIndexCurrentCellDensity)   = dxInv*(FL(0) - FR(0));
      V(vIndexCurrentCellDensity+1) = dxInv*(FL(1) - FR(1));
      V(vIndexCurrentCellDensity+2) = dxInv*(FL(2) - FR(2));

      if(J){
	CellJacobianFunctor(smPt, bcCellJacFactors);
      }
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

  void allocateGhostValues()
  {
    /*
      In 1d, the ghost values at left and right are stored in arrays.
      The size of this array depends on the stencil.
      For now, we has as many rows as BD cells so that for each
      BD cell, the corresponding ghosts can be easily found.
      This could be changed since in 1d we literally just need two arrays only.

      stencil = 3:
      ghost array has 3 entries for the 3 dofs in the ghost cell.
      For example at left boundary we have:
      ----------------||
      |	 0,  1,     2 ||
      | rho, rho*u, E ||
      |		      ||
      ----------------||

      stencil = 5:
      ghost array has 6 entries such that at left boundary we have:
      note the order of the indexing
      ----------------|---------------||
      |	 3,  4,     5 |	 0,  1,     2 ||
      | rho, rho*u, E | rho, rho*u, E ||
      |		      |		      ||
      ----------------|---------------||

      stencil = 7:
      ghost array has 9 entries such that at left boundary we have:
      note the order of the indexing
      ----------------|---------------|---------------||
      |	 6,  7,     8 |	 3,  4,     5 |	 0,  1,     2 ||
      | rho, rho*u, E | rho, rho*u, E | rho, rho*u, E ||
      |		      |		      |		      ||
      ----------------|---------------|---------------||
     */

    const auto stencilSize    = reconstructionTypeToStencilSize(m_recEn);
    const auto numGhostValues = numDofPerCell*((stencilSize-1)/2);
    const index_t s1 = m_meshObj.numCellsBd();
    ::pressiodemoapps::resize(m_ghostLeft,  s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight, s1, numGhostValues);
  }

protected:
  scalar_type m_gamma = static_cast<scalar_type>(1.4);

  const MeshType & m_meshObj;
  ::pressiodemoapps::Euler1d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  ::pressiodemoapps::InviscidFluxScheme m_fluxEn;

  // note that dof refers to the degress of freedom,
  // which is NOT same as grid points. for this problem,
  // the dof = numDofPerCell * number_of_unknown_grid_points
  // SampleMesh_ identifies the velocity/residual locations
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  mutable stencil_container_type m_stencilVals = {};
  mutable ghost_container_type   m_ghostLeft = {};
  mutable ghost_container_type   m_ghostRight = {};
};

}}
#endif
