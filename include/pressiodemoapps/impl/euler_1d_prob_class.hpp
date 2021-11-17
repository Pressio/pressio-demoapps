
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
#include "euler_flux_functor.hpp"
#include "functor_cell_jacobian.hpp"
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
  using ghost_container_type      = Eigen::Matrix<scalar_type,-1,-1, Eigen::RowMajor>;
  using stencil_container_type    = Eigen::Matrix<scalar_type,-1,1>;
  using flux_type	          = Eigen::Matrix<scalar_type,-1,1>;
  using edge_rec_type	          = Eigen::Matrix<scalar_type,-1,1>;
  using flux_jac_type             = Eigen::Matrix<scalar_type, numDofPerCell, numDofPerCell>;
  using reconstruction_gradient_t = Eigen::Matrix<scalar_type, -1, -1>;

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

    using Tr = Eigen::Triplet<scalar_type>;
    std::vector<Tr> trList;
    initializeJacobianForNearBoundaryCells(trList);
    initializeJacobianForInnerCells(trList);
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

    // reconstructions values
    edge_rec_type uMinusHalfNeg(numDofPerCell);
    edge_rec_type uMinusHalfPos(numDofPerCell);
    edge_rec_type uPlusHalfNeg(numDofPerCell);
    edge_rec_type uPlusHalfPos(numDofPerCell);
    // fluxes
    flux_type fluxL(numDofPerCell);
    flux_type fluxR(numDofPerCell);
    // flux jacobians
    flux_jac_type fluxJacLNeg(numDofPerCell, numDofPerCell);
    flux_jac_type fluxJacLPos(numDofPerCell, numDofPerCell);
    flux_jac_type fluxJacRNeg(numDofPerCell, numDofPerCell);
    flux_jac_type fluxJacRPos(numDofPerCell, numDofPerCell);

    fillGhostsIfNeeded(U, currentTime);

    if (J){
      ::pressiodemoapps::set_zero(*J);
    }

    int nonZerosCountBeforeComputing = 0;
    if (J){
      nonZerosCountBeforeComputing = J->nonZeros();
      ::pressiodemoapps::set_zero(*J);
    }

    velocityAndJacInnerCellsImpl(U, currentTime, V, J,
				 fluxL, fluxR,
				 fluxJacLNeg, fluxJacLPos,
				 fluxJacRNeg, fluxJacRPos,
				 uMinusHalfNeg, uMinusHalfPos,
				 uPlusHalfNeg,  uPlusHalfPos);

    velocityAndJacNearBDCellsImpl(U, currentTime, V, J,
				  fluxL, fluxR,
				  fluxJacLNeg, fluxJacLPos,
				  fluxJacRNeg, fluxJacRPos,
				  uMinusHalfNeg, uMinusHalfPos,
				  uPlusHalfNeg,  uPlusHalfPos);
  }

private:
  template<class Tr>
  void initializeJacobianForInnerCells(std::vector<Tr> & trList)
  {
    const scalar_type zero = 0;
    const auto & graph = m_meshObj.graph();
    const auto & targetGraphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (std::size_t it=0; it<targetGraphRows.size(); ++it)
      {
	const auto smPt = targetGraphRows[it];
	const auto jacRowOfCurrCellRho = smPt*numDofPerCell;

	// entries wrt current cell's dofs
	const auto jacColOfCurrCellRho = graph(smPt, 0)*numDofPerCell;
	for (int k=0; k<numDofPerCell; ++k){
	  for (int j=0; j<numDofPerCell; ++j){
	    trList.push_back( Tr(jacRowOfCurrCellRho+k, jacColOfCurrCellRho+j, zero) );
	  }
	}

	// wrt neighbors
	const auto numNeighbors =
	  (m_recEn == InviscidFluxReconstruction::FirstOrder) ? 2
	    : (m_recEn == InviscidFluxReconstruction::Weno3) ? 4 : 6;

	for (int i=1; i<=numNeighbors; ++i){
	  const auto ci = graph(smPt, i)*numDofPerCell;
	  for (int k=0; k<numDofPerCell; ++k){
	    for (int j=0; j<numDofPerCell; ++j){
	      trList.push_back( Tr(jacRowOfCurrCellRho+k, ci+j, zero) );
	    }
	  }
	}
      }
  }

  template<class Tr>
  void initializeJacobianForNearBoundaryCells(std::vector<Tr> & trList)
  {
    const scalar_type zero = 0;
    const auto & graph = m_meshObj.graph();
    const auto & targetGraphRows = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<targetGraphRows.size(); ++it)
      {
	const auto smPt = targetGraphRows[it];
	const auto jacRowOfCurrCellRho = smPt*numDofPerCell;
	const auto jacColOfCurrCellRho = graph(smPt, 0)*numDofPerCell;

	// wrt current cell's dofs
	for (int k=0; k<numDofPerCell; ++k){
	  for (int j=0; j<numDofPerCell; ++j){
	    trList.push_back( Tr(jacRowOfCurrCellRho+k, jacColOfCurrCellRho+j, zero) );
	  }
	}

	// wrt neighbors
	// for near-bd, we only do first-order Jacobian for now
	const auto L0 = graph(smPt, 1);
	if (L0 != -1){
	  const auto ciL0 = L0*numDofPerCell;
	  for (int k=0; k<numDofPerCell; ++k){
	    for (int j=0; j<numDofPerCell; ++j){
	      trList.push_back( Tr(jacRowOfCurrCellRho+k, ciL0+j, zero) );
	    }
	  }
	}

	const auto R0 = graph(smPt, 2);
	if (R0 != -1){
	  const auto ciR0 = R0*numDofPerCell;
	  for (int k=0; k<numDofPerCell; ++k){
	    for (int j=0; j<numDofPerCell; ++j){
	      trList.push_back( Tr(jacRowOfCurrCellRho+k, ciR0+j, zero) );
	    }
	  }
	}
      }
  }

  template<class U_t>
  void fillGhostsIfNeeded(const U_t & U,
			  scalar_type /*currTime*/) const
  {
    namespace pda = ::pressiodemoapps;

    const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_recEn);
    // only need ghosts for specific problems
    if (m_probEn == pda::Euler1d::Sod or m_probEn == pda::Euler1d::Lax)
    {
      using ghost_filler_t  = pda::impl::Ghost1dNeumannFiller<
	numDofPerCell, U_t, MeshType, ghost_container_type>;
      ghost_filler_t ghF(stencilSizeNeeded, U, m_meshObj, m_ghostLeft, m_ghostRight);

      const auto & rowsBd = m_meshObj.graphRowsOfCellsNearBd();
      for (int it=0; it<rowsBd.size(); ++it){
	ghF(rowsBd[it], it);
      }
    }
  }

  template<class U_t, class V_t>
  void velocityAndJacInnerCellsImpl(const U_t & U,
				    const scalar_type currentTime,
				    V_t & V,
				    jacobian_type * J,
				    flux_type & fluxL,
				    flux_type & fluxR,
				    flux_jac_type & fluxJacLNeg,
				    flux_jac_type & fluxJacLPos,
				    flux_jac_type & fluxJacRNeg,
				    flux_jac_type & fluxJacRPos,
				    edge_rec_type & uMinusHalfNeg,
				    edge_rec_type & uMinusHalfPos,
				    edge_rec_type & uPlusHalfNeg,
				    edge_rec_type & uPlusHalfPos) const
  {
    // for inner cells, both velocity and Jacobian
    // are computed according to the order selected by the user

    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;

    // allocate gradients of reconstructed states
    // the size depends on the scheme selected
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    reconstruction_gradient_t gradLNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradLPos(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRPos(numDofPerCell, stencilSize-1);

    if (J){

      // using functor_type =
      // 	Velocity<
      // 	  Jacobian<
      // 	    CellFluxesAndGradients<
      // 	      UsingReconstructionFromStateWithGrad<
      // 		dimensionality, numDofPerCell, edge_rec_type, U_t, MeshType, reconstruction_gradient_t>,
      // 	      flux_type, flux_jac_type>,
      // 	    jacobian_type>,
      // 	  U_t
      // 	>;

      // functor_type F(U, *J,
      // 		     fluxL, fluxR, fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
      // 		     xAxis, m_recEn, m_meshObj, U,
      // 		     uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
      // 		     gradLNeg, gradLPos, gradRNeg, gradRPos);

      using reconstruct_functor_t = ::pressiodemoapps::impl::ReconstructorFromState<
	dimensionality, numDofPerCell, edge_rec_type, U_t,
	MeshType, reconstruction_gradient_t>;
      reconstruct_functor_t Reconstructor(xAxis, m_recEn, U, m_meshObj,
					  uMinusHalfNeg, uMinusHalfPos,
					  uPlusHalfNeg,  uPlusHalfPos,
					  gradLNeg, gradLPos,
					  gradRNeg, gradRPos);

      using flux_fnct_t = pda::ee::impl::FluxFunctorWithGradients<
	numDofPerCell, scalar_type, edge_rec_type, flux_type, flux_jac_type>;
      flux_fnct_t FluxFunctor(m_fluxEn, m_gamma,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos,
			      fluxL, fluxR,
			      fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos);

      using jac_fnct_t = ::pressiodemoapps::impl::InnerCellJacobianFunctor<
	dimensionality, numDofPerCell, scalar_type, jacobian_type,
	flux_jac_type, MeshType, reconstruction_gradient_t>;
      jac_fnct_t CellJacobianFunctor(m_recEn, *J, m_meshObj,
				     fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
				     gradLNeg, gradLPos, gradRNeg, gradRPos);

      innerCellsLoopImpl(U, currentTime, V, J, fluxL, fluxR,
			 Reconstructor, FluxFunctor, CellJacobianFunctor);
    }
    else{

      using reconstruct_functor_t = ::pressiodemoapps::impl::ReconstructorFromState<
	dimensionality, numDofPerCell, edge_rec_type, U_t, MeshType>;
      reconstruct_functor_t Reconstructor(xAxis, m_recEn, U, m_meshObj,
					  uMinusHalfNeg, uMinusHalfPos,
					  uPlusHalfNeg,  uPlusHalfPos);

      using flux_fnct_t = pda::ee::impl::FluxFunctorWithoutGradients<
	numDofPerCell, scalar_type, edge_rec_type, flux_type>;
      flux_fnct_t FluxFunctor(m_fluxEn, m_gamma,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos,
			      fluxL, fluxR);

      innerCellsLoopImpl(U, currentTime, V, J,fluxL, fluxR,
			 Reconstructor, FluxFunctor,
			 // no-op for cell jacobian since Jacobian is disabled
			 [](index_t /*unused*/){});
    }

  }

  template<
    class U_t, class V_t,
    class ReconstructFunctor,
    class FluxFunctor,
    class CellJacobianFunctor
    >
  void innerCellsLoopImpl(const U_t & U,
			  const scalar_type currentTime,
			  V_t & V,
			  jacobian_type * J,
			  const flux_type & fluxL,
			  const flux_type & fluxR,
			  ReconstructFunctor & reconstructor,
			  FluxFunctor & fluxesWithOptionalGradsFunctor,
			  CellJacobianFunctor && cellJacobianFunctor) const
  {
    const auto dxInv = m_meshObj.dxInv();
    const auto & graphRows = m_meshObj.graphRowsOfCellsAwayFromBd();

    for (std::size_t it=0; it<graphRows.size(); ++it)
      {
	const auto smPt = graphRows[it];
	reconstructor(smPt);
	fluxesWithOptionalGradsFunctor();
	const auto vIndexCurrentCellDensity = smPt*numDofPerCell;
	V(vIndexCurrentCellDensity)   = dxInv*(fluxL(0) - fluxR(0));
	V(vIndexCurrentCellDensity+1) = dxInv*(fluxL(1) - fluxR(1));
	V(vIndexCurrentCellDensity+2) = dxInv*(fluxL(2) - fluxR(2));

	if(J){
	  cellJacobianFunctor(smPt);
	}
      }
  }

  template<class U_t, class V_t>
  void velocityAndJacNearBDCellsImpl(const U_t & U,
				     const scalar_type currentTime,
				     V_t & V,
				     jacobian_type * J,
				     flux_type & fluxL,
				     flux_type & fluxR,
				     flux_jac_type & fluxJacLNeg,
				     flux_jac_type & fluxJacLPos,
				     flux_jac_type & fluxJacRNeg,
				     flux_jac_type & fluxJacRPos,
				     edge_rec_type & uMinusHalfNeg,
				     edge_rec_type & uMinusHalfPos,
				     edge_rec_type & uPlusHalfNeg,
				     edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;

    // for near-BD cells, scheme for velocity is whatever user decides,
    // but Jacobian is limited to firstorder

    // allocate gradients of reconstructed states (MUST BE FIRST-ORDER)
    const auto jacRecEnum = pda::InviscidFluxReconstruction::FirstOrder;
    const auto stencilSizeForJac = reconstructionTypeToStencilSize(jacRecEnum);
    reconstruction_gradient_t gradLNeg(numDofPerCell, stencilSizeForJac-1);
    reconstruction_gradient_t gradLPos(numDofPerCell, stencilSizeForJac-1);
    reconstruction_gradient_t gradRNeg(numDofPerCell, stencilSizeForJac-1);
    reconstruction_gradient_t gradRPos(numDofPerCell, stencilSizeForJac-1);

    if (m_recEn == InviscidFluxReconstruction::FirstOrder)
    {
      // if here, then the scheme for velocity matches
      // the one for Jacobian so we can use same functors

      using stencil_filler_t  = pda::impl::StencilFiller<
	dimensionality, numDofPerCell, stencil_container_type,
	U_t, MeshType, ghost_container_type>;
      stencil_filler_t FillStencilFunctor(reconstructionTypeToStencilSize(m_recEn), U,
					  m_meshObj, m_ghostLeft, m_ghostRight,
					  m_stencilVals);

      using rec_functor_t = pda::impl::ReconstructorFromStencilThreeDofPerCell<
	edge_rec_type, stencil_container_type>;
      rec_functor_t ReconstructFunctor(m_recEn, m_stencilVals,
				       uMinusHalfNeg, uMinusHalfPos,
				       uPlusHalfNeg,  uPlusHalfPos);

      using flux_fnct_t = pda::ee::impl::FluxFunctorWithGradients<
	numDofPerCell, scalar_type, edge_rec_type, flux_type, flux_jac_type>;
      flux_fnct_t FluxFunctor(m_fluxEn, m_gamma,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos,
			      fluxL, fluxR,
			      fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos);

      using jac_functor_t = pda::impl::FirstOrderBdCellJacobianFunctor<
	dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, MeshType>;
      jac_functor_t CellJacobianFunctor(*J, m_meshObj, fluxJacLNeg, fluxJacLPos,
					fluxJacRNeg, fluxJacRPos);

      nearBDCellsLoopImpl(U, currentTime, V, J, fluxL, fluxR,
			  FillStencilFunctor,
			  // noop lambda because velocity and jacobian use same filler
			  [](index_t /**/, index_t /**/){},
			  ReconstructFunctor,
			  // noop because velocity and jacobian use same reconstuctor
			  [](){},
			  FluxFunctor,
			  CellJacobianFunctor);
    }

    else if (m_recEn == InviscidFluxReconstruction::Weno3 or
	     m_recEn == InviscidFluxReconstruction::Weno5)
    {
      // if here, then for the fluxes we need to meet the
      // scheme desired by user, while for jacobian we keep firstorder

      using stencil_filler_t  = pda::impl::StencilFiller<
	dimensionality, numDofPerCell, stencil_container_type,
	U_t, MeshType, ghost_container_type>;

      using rec_functor_t = pda::impl::ReconstructorFromStencilThreeDofPerCell<
	edge_rec_type, stencil_container_type>;

      // instantiate functors for computing velocity
      stencil_filler_t FillStencilFunctor(reconstructionTypeToStencilSize(m_recEn), U,
					  m_meshObj, m_ghostLeft, m_ghostRight,
					  m_stencilVals);
      rec_functor_t ReconstructFunctor(m_recEn, m_stencilVals,
				       uMinusHalfNeg, uMinusHalfPos,
				       uPlusHalfNeg,  uPlusHalfPos);

      // for jacobian i need auxiliary stuff, because Jacobian is first order
      // which is not same scheme wanted for velocity
      edge_rec_type uMinusHalfNegForJ(numDofPerCell);
      edge_rec_type uMinusHalfPosForJ(numDofPerCell);
      edge_rec_type uPlusHalfNegForJ(numDofPerCell);
      edge_rec_type uPlusHalfPosForJ(numDofPerCell);

      const auto cellJacOrdEn = pda::InviscidFluxReconstruction::FirstOrder;
      const auto stencilSizeForJ = reconstructionTypeToStencilSize(cellJacOrdEn);
      stencil_container_type stencilValsForJ(numDofPerCell*stencilSizeForJ);
      stencil_filler_t FillStencilFunctorForJ(stencilSizeForJ, U, m_meshObj,
					      m_ghostLeft, m_ghostRight,
					      stencilValsForJ);

      rec_functor_t ReconstructFunctorForJ(cellJacOrdEn, stencilValsForJ,
					   uMinusHalfNegForJ, uMinusHalfPosForJ,
					   uPlusHalfNegForJ,  uPlusHalfPosForJ);

      using flux_fnct_t = pda::ee::impl::FluxFunctorGradUsingDifferentRecontructions<
	numDofPerCell, scalar_type, edge_rec_type, flux_type, flux_jac_type>;
      flux_fnct_t FluxFunctor(m_fluxEn, m_gamma,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg, uPlusHalfPos,
			      uMinusHalfNegForJ, uMinusHalfPosForJ,
			      uPlusHalfNegForJ, uPlusHalfPosForJ,
			      fluxL, fluxR,
			      fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos);

      using jac_functor_t = pda::impl::FirstOrderBdCellJacobianFunctor<
	dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, MeshType>;
      jac_functor_t CellJacobianFunctor(*J, m_meshObj, fluxJacLNeg, fluxJacLPos,
					fluxJacRNeg, fluxJacRPos);

      nearBDCellsLoopImpl(U, currentTime, V, J, fluxL, fluxR,
			  FillStencilFunctor, FillStencilFunctorForJ,
			  ReconstructFunctor, ReconstructFunctorForJ,
			  FluxFunctor, CellJacobianFunctor);
    }
  }

  template<
    class U_t, class V_t,
    class StencilFillFunctorForV,
    class StencilFillFunctorForJ,
    class ReconstructFunctorForV,
    class ReconstructFunctorForJ,
    class FluxFunctor,
    class CellJacobianFunctor
    >
  void nearBDCellsLoopImpl(const U_t & U,
			   const scalar_type currentTime,
			   V_t & V,
			   jacobian_type * J,
			   flux_type & fluxL,
			   flux_type & fluxR,
			   StencilFillFunctorForV & stencilFillerForVelo,
			   StencilFillFunctorForJ && stencilFillerForJac,
			   ReconstructFunctorForV & reconstructorForVelo,
			   ReconstructFunctorForJ && reconstructorForJac,
			   FluxFunctor & fluxFunctor,
			   CellJacobianFunctor && cellJacobianFunctor) const
  {
    namespace pda = ::pressiodemoapps;

    std::array<scalar_type, numDofPerCell> bdCellJacFactors;
    bdCellJacFactors.fill(static_cast<scalar_type>(1));

    const auto dxInv = m_meshObj.dxInv();
    const auto & graphRows = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];

      stencilFillerForVelo(smPt, it);
      reconstructorForVelo();
      stencilFillerForJac(smPt, it);
      reconstructorForJac();
      fluxFunctor();
      const auto vIndexCurrentCellDensity = smPt*numDofPerCell;
      V(vIndexCurrentCellDensity)   = dxInv*(fluxL(0) - fluxR(0));
      V(vIndexCurrentCellDensity+1) = dxInv*(fluxL(1) - fluxR(1));
      V(vIndexCurrentCellDensity+2) = dxInv*(fluxL(2) - fluxR(2));

      if(J){
	cellJacobianFunctor(smPt, bdCellJacFactors);
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

} // end namespace implee1d
} // end namespace pressiodemoapps
#endif






  //   // reconstruct functor for face fluxes
  //   // here we need to use whatever order (m_recEn) user decides
  //   using rec_fnct_t = pda::impl::ReconstructorFromState<
  //     dimensionality, edge_rec_type, U_t, MeshType>;
  //   rec_fnct_t Reconstructor(m_recEn, U, m_meshObj,
  // 			     uMinusHalfNeg, uMinusHalfPos,
  // 			     uPlusHalfNeg,  uPlusHalfPos);

  //   // cell jacobian functor
  //   // currently, REGARDLESS of the reconstruction scheme,
  //   // we have only first-order Jacobian so we need a first-order reconstructor
  //   // for the Jacobian
  //   edge_rec_type uMinusHalfNegForJ(numDofPerCell);
  //   edge_rec_type uMinusHalfPosForJ(numDofPerCell);
  //   edge_rec_type uPlusHalfNegForJ(numDofPerCell);
  //   edge_rec_type uPlusHalfPosForJ(numDofPerCell);
  //   rec_fnct_t ReconstructorForJ(pda::InviscidFluxReconstruction::FirstOrder,
  // 				 U, m_meshObj,
  // 				 uMinusHalfNegForJ, uMinusHalfPosForJ,
  // 				 uPlusHalfNegForJ,  uPlusHalfPosForJ);

  //   using jac_fnct_t = pda::impl::FirstOrderInnerCellJacobianFunctor<
  //     dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, MeshType>;
  //   jac_fnct_t CellJacobianFunctor(*J, m_meshObj, JLneg, JLpos, JRneg, JRpos);

  //   // loop
  //   const auto dxInv = m_meshObj.dxInv();
  //   const auto & targetGraphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
  //   for (int it=0; it<targetGraphRows.size(); ++it)
  //     {
  // 	const auto smPt = targetGraphRows[it];
  // 	Reconstructor.template operator()<numDofPerCell>(smPt);

  // 	if(J){
  // 	  // note that REGARDLESS of the reconstruction scheme,
  // 	  // we currently only have only first-order Jacobian so we need
  // 	  // to run the reconstructor for the Jacobian
  // 	  // which will ensure that uMinusNegForJ, etc have the right values
  // 	  ReconstructorForJ.template operator()<numDofPerCell>(smPt);
  // 	}

  // 	switch(m_fluxEn)
  // 	  {
  // 	  case pda::InviscidFluxScheme::Rusanov:
  // 	    ee::impl::eeRusanovFluxThreeDof(FL, uMinusHalfNeg, uMinusHalfPos, m_gamma);
  // 	    ee::impl::eeRusanovFluxThreeDof(FR, uPlusHalfNeg,  uPlusHalfPos,  m_gamma);

  // 	    if(J){
  // 	      ee::impl::eeRusanovFluxJacobianThreeDof(JLneg, JLpos, uMinusHalfNegForJ,
  // 						      uMinusHalfPosForJ, m_gamma);
  // 	      ee::impl::eeRusanovFluxJacobianThreeDof(JRneg, JRpos, uPlusHalfNegForJ,
  // 						      uPlusHalfPosForJ, m_gamma);
  // 	    }
  // 	    break;
  // 	  }

  // 	const auto vIndex = smPt*numDofPerCell;
  // 	V(vIndex)   = dxInv*(FL(0) - FR(0));
  // 	V(vIndex+1) = dxInv*(FL(1) - FR(1));
  // 	V(vIndex+2) = dxInv*(FL(2) - FR(2));

  // 	if(J){
  // 	  CellJacobianFunctor(smPt);
  // 	}
  //     }
  // }


  // template<class U_t, class V_t>
  // void velocityAndJacNearBDCellsImpl(const U_t & U,
  // 				     const scalar_type currentTime,
  // 				     V_t & V,
  // 				     jacobian_type * J,
  // 				     flux_type & FL,
  // 				     flux_type & FR,
  // 				     flux_jac_type & JLneg,
  // 				     flux_jac_type & JLpos,
  // 				     flux_jac_type & JRneg,
  // 				     flux_jac_type & JRpos,
  // 				     edge_rec_type & uMinusHalfNeg,
  // 				     edge_rec_type & uMinusHalfPos,
  // 				     edge_rec_type & uPlusHalfNeg,
  // 				     edge_rec_type & uPlusHalfPos) const
  // {
  //   namespace pda = ::pressiodemoapps;

  //   using stencil_filler_t  = pda::impl::StencilFiller<
  //     dimensionality, numDofPerCell, stencil_container_type,
  //     U_t, MeshType, ghost_container_type>;

  //   using rec_functor_t = pda::impl::ReconstructorFromStencilThreeDofPerCell<
  //     edge_rec_type, stencil_container_type>;

  //   using jac_functor_t = pda::impl::FirstOrderBdCellJacobianFunctor<
  //     dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, MeshType>;

  //   // -----------------------------------------------
  //   // functors needed to compute face fluxes
  //   // NOTE: for face fluxes, we must be consistent with whatever scheme
  //   // was specified by user, so we need to use this enum: m_recEn
  //   // -----------------------------------------------
  //   stencil_filler_t FillStencilValuesFunctor(reconstructionTypeToStencilSize(m_recEn), U,
  // 					      m_meshObj, m_ghostLeft, m_ghostRight,
  // 					      m_stencilVals);
  //   rec_functor_t FaceValuesReconstructFunctor(m_recEn, m_stencilVals,
  // 					       uMinusHalfNeg, uMinusHalfPos,
  // 					       uPlusHalfNeg,  uPlusHalfPos);

  //   // ----------------------------------------
  //   // functors needed to compute cell jacobian
  //   // we have only first-order Jacobian for now, so we need
  //   // dedicate functors because we cannot use those above.
  //   // Once we have cell Jacobians of various order, we can change this.
  //   // ----------------------------------------
  //   const auto cellJacOrdEn = pda::InviscidFluxReconstruction::FirstOrder;

  //   stencil_container_type stencilValsForJ = {};
  //   const auto stencilSizeForJ = reconstructionTypeToStencilSize(cellJacOrdEn);
  //   pda::resize(stencilValsForJ, numDofPerCell*stencilSizeForJ);
  //   stencil_filler_t FillStencilValuesFunctorForJ(stencilSizeForJ, U, m_meshObj,
  // 						  m_ghostLeft, m_ghostRight, stencilValsForJ);

  //   edge_rec_type uMinusHalfNegForJ(numDofPerCell);
  //   edge_rec_type uMinusHalfPosForJ(numDofPerCell);
  //   edge_rec_type uPlusHalfNegForJ(numDofPerCell);
  //   edge_rec_type uPlusHalfPosForJ(numDofPerCell);
  //   rec_functor_t FaceValuesReconstructFunctorForJ(cellJacOrdEn, stencilValsForJ,
  // 						   uMinusHalfNegForJ, uMinusHalfPosForJ,
  // 						   uPlusHalfNegForJ,  uPlusHalfPosForJ);
  //   jac_functor_t CellJacobianFunctor(*J, m_meshObj, JLneg, JLpos, JRneg, JRpos);

  //   std::array<scalar_type, numDofPerCell> bcCellJacFactors;
  //   bcCellJacFactors.fill(static_cast<scalar_type>(1));

  //   // ----------------------------------------
  //   // loop over cells
  //   // ----------------------------------------
  //   const auto dxInv = m_meshObj.dxInv();
  //   const auto & targetGraphRows = m_meshObj.graphRowsOfCellsNearBd();
  //   for (std::size_t it=0; it<targetGraphRows.size(); ++it)
  //   {
  //     const auto smPt = targetGraphRows[it];

  //     FillStencilValuesFunctor(smPt, it);
  //     FaceValuesReconstructFunctor();

  //     if(J){
  // 	// note that REGARDLESS of the reconstruction scheme,
  // 	// we currently only have only first-order Jacobian so we need
  // 	// to run the reconstructor for the Jacobian
  // 	// which will ensure that uMinusNegForJ, etc have the right values
  // 	FillStencilValuesFunctorForJ(smPt, it);
  // 	FaceValuesReconstructFunctorForJ();
  //     }

  //     switch(m_fluxEn)
  //     {
  //     case pda::InviscidFluxScheme::Rusanov:
  // 	ee::impl::eeRusanovFluxThreeDof(FL, uMinusHalfNeg, uMinusHalfPos, m_gamma);
  // 	ee::impl::eeRusanovFluxThreeDof(FR, uPlusHalfNeg,  uPlusHalfPos,  m_gamma);

  // 	if (J){
  // 	  ee::impl::eeRusanovFluxJacobianThreeDof(JLneg, JLpos, uMinusHalfNegForJ,
  // 						  uMinusHalfPosForJ, m_gamma);
  // 	  ee::impl::eeRusanovFluxJacobianThreeDof(JRneg, JRpos, uPlusHalfNegForJ,
  // 						  uPlusHalfPosForJ, m_gamma);
  // 	}
  // 	break;
  //     }

  //     const auto vIndexCurrentCellDensity = smPt*numDofPerCell;
  //     V(vIndexCurrentCellDensity)   = dxInv*(FL(0) - FR(0));
  //     V(vIndexCurrentCellDensity+1) = dxInv*(FL(1) - FR(1));
  //     V(vIndexCurrentCellDensity+2) = dxInv*(FL(2) - FR(2));

  //     if(J){
  // 	CellJacobianFunctor(smPt, bcCellJacFactors);
  //     }
  //   }
  // }
