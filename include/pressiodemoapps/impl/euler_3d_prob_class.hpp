
#ifndef PRESSIODEMOAPPS_EULER3D_HPP_
#define PRESSIODEMOAPPS_EULER3D_HPP_

#include "euler_rankine_hugoniot.hpp"
#include "euler_rusanov_flux_values_function.hpp"
#include "euler_rusanov_flux_jacobian_function.hpp"
#include "euler_3d_initial_condition.hpp"
#include "functor_ghost_fill_neumann.hpp"
#include "euler_3d_ghost_filler.hpp"
#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "euler_flux_mixin.hpp"
#include "mixin_directional_flux_balance.hpp"
#include "mixin_directional_flux_balance_jacobian.hpp"

namespace pressiodemoapps{ namespace ee{ namespace impl{

template<class MeshType>
class EigenEuler3dApp
{

public:
  using index_t	      = typename MeshType::index_t;
  using scalar_type   = typename MeshType::scalar_t;
  using state_type    = Eigen::Matrix<scalar_type,-1,1>;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, index_t>;

  static constexpr int	   dimensionality{3};
  static constexpr index_t numDofPerCell{5};

private:
  using ghost_container_type      = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using stencil_container_type    = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using flux_type	          = Eigen::Matrix<scalar_type, numDofPerCell,  1>;
  using edge_rec_type	          = Eigen::Matrix<scalar_type, numDofPerCell,  1>;
  using flux_jac_type             = Eigen::Matrix<scalar_type, numDofPerCell, numDofPerCell>;
  using reconstruction_gradient_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;

public:
  EigenEuler3dApp(const MeshType & meshObj,
		  ::pressiodemoapps::Euler3d probEnum,
		  ::pressiodemoapps::InviscidFluxReconstruction recEnum,
		  ::pressiodemoapps::InviscidFluxScheme fluxEnum)
    : m_meshObj(meshObj), m_probEn(probEnum),
      m_recEn(recEnum), m_fluxEn(fluxEnum)
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

protected:
  void initializeJacobian(jacobian_type & J)
  {
    J.resize(m_numDofSampleMesh, m_numDofStencilMesh);

    using Tr = Eigen::Triplet<scalar_type>;
    std::vector<Tr> trList;
    initializeJacobianForNearBoundaryCells(trList);
    initializeJacobianForInnerCells(trList);
    J.setFromTriplets(trList.begin(), trList.end());
    // compress to make it Csr
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
    edge_rec_type uMinusHalfNeg, uMinusHalfPos;
    edge_rec_type uPlusHalfNeg,  uPlusHalfPos;
    // fluxes
    flux_type fluxL, fluxR; //left, right x
    flux_type fluxB, fluxF; //back, front y
    flux_type fluxD, fluxU; //down, up z

    // flux jacobians
    flux_jac_type fluxJacLNeg, fluxJacLPos;
    flux_jac_type fluxJacRNeg, fluxJacRPos;
    flux_jac_type fluxJacBNeg, fluxJacBPos;
    flux_jac_type fluxJacFNeg, fluxJacFPos;
    flux_jac_type fluxJacDNeg, fluxJacDPos;
    flux_jac_type fluxJacUNeg, fluxJacUPos;

    fillGhosts(U);

    V.setZero();

    int nonZerosCountBeforeComputing = 0;
    if (J){
      nonZerosCountBeforeComputing = J->nonZeros();
      ::pressiodemoapps::set_zero(*J);

      // near boundary I have be careful because
      // the jacobian can only be first order for now
      if (m_recEn == InviscidFluxReconstruction::FirstOrder){
	velocityAndJacNearBDCellsImplFirstOrder(U, currentTime, V, *J,
						fluxL, fluxR, fluxB, fluxF, fluxD, fluxU,
						fluxJacLNeg, fluxJacLPos,
						fluxJacRNeg, fluxJacRPos,
						fluxJacBNeg, fluxJacBPos,
						fluxJacFNeg, fluxJacFPos,
						fluxJacDNeg, fluxJacDPos,
						fluxJacUNeg, fluxJacUPos,
						uMinusHalfNeg, uMinusHalfPos,
						uPlusHalfNeg,  uPlusHalfPos);
      }

      else{
	velocityAndJacNearBDCellsImplDifferentScheme(U, currentTime, V, *J,
						     fluxL, fluxR, fluxB, fluxF, fluxD, fluxU,
						     fluxJacLNeg, fluxJacLPos,
						     fluxJacRNeg, fluxJacRPos,
						     fluxJacBNeg, fluxJacBPos,
						     fluxJacFNeg, fluxJacFPos,
						     fluxJacDNeg, fluxJacDPos,
						     fluxJacUNeg, fluxJacUPos,
						     uMinusHalfNeg, uMinusHalfPos,
						     uPlusHalfNeg,  uPlusHalfPos);
      }

      velocityAndJacInnerCellsImpl(U, currentTime, V, *J,
				   fluxL, fluxR, fluxB, fluxF, fluxD, fluxU,
				   fluxJacLNeg, fluxJacLPos,
				   fluxJacRNeg, fluxJacRPos,
				   fluxJacBNeg, fluxJacBPos,
				   fluxJacFNeg, fluxJacFPos,
				   fluxJacDNeg, fluxJacDPos,
				   fluxJacUNeg, fluxJacUPos,
				   uMinusHalfNeg, uMinusHalfPos,
				   uPlusHalfNeg,  uPlusHalfPos);

      assert(J->nonZeros() == nonZerosCountBeforeComputing);
    }

    else{
      velocityOnlyNearBdCellsImpl(U, currentTime, V,
				  fluxL, fluxR, fluxB, fluxF, fluxD, fluxU,
				  uMinusHalfNeg, uMinusHalfPos,
				  uPlusHalfNeg,  uPlusHalfPos);

      velocityOnlyInnerCellsImpl(U, currentTime, V,
				 fluxL, fluxR, fluxB, fluxF, fluxD, fluxU,
				 uMinusHalfNeg, uMinusHalfPos,
				 uPlusHalfNeg,  uPlusHalfPos);
    }
  }

private:
  template<class Tr>
  void initializeJacobianForInnerCells(std::vector<Tr> & trList)
  {
    // for inner cells, the Jacobian is of the same scheme
    // wanted by the user, no special treatment is needed

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

	// wrt neighbors: this depends on the advection scheme
	const auto numNeighbors =
	  (m_recEn == InviscidFluxReconstruction::FirstOrder) ? 6
	    : (m_recEn == InviscidFluxReconstruction::Weno3) ? 12 : 18;

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
	for (int i=1; i<=6; ++i){
	  const auto nID = graph(smPt, i);
	  if (nID != -1){
	    const auto ci = nID*numDofPerCell;
	    for (int k=0; k<numDofPerCell; ++k){
	      for (int j=0; j<numDofPerCell; ++j){
		trList.push_back( Tr(jacRowOfCurrCellRho+k, ci+j, zero) );
	      }
	    }
	  }
	}
      }
  }

  template<class U_t>
  void fillGhosts(const U_t & U) const
  {
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    if (m_probEn == ::pressiodemoapps::Euler3d::SedovSymmetry)
    {
      using ghost_filler_t  = ::pressiodemoapps::ee::impl::Ghost3dSedov<U_t, MeshType, ghost_container_type>;
      ghost_filler_t ghF(stencilSize, U, m_meshObj,
			 m_ghostLeft, m_ghostRight,
			 m_ghostBack, m_ghostFront,
			 m_ghostDown, m_ghostUp);

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

  template<class U_t, class V_t>
  void velocityAndJacInnerCellsImpl(const U_t & U,
				    const scalar_type currentTime,
				    V_t & V,
				    jacobian_type & J,
				    flux_type & fluxL,
				    flux_type & fluxR,
				    flux_type & fluxB,
				    flux_type & fluxF,
				    flux_type & fluxD,
				    flux_type & fluxU,
				    flux_jac_type & fluxJacLNeg,
				    flux_jac_type & fluxJacLPos,
				    flux_jac_type & fluxJacRNeg,
				    flux_jac_type & fluxJacRPos,
				    flux_jac_type & fluxJacBNeg,
				    flux_jac_type & fluxJacBPos,
				    flux_jac_type & fluxJacFNeg,
				    flux_jac_type & fluxJacFPos,
				    flux_jac_type & fluxJacDNeg,
				    flux_jac_type & fluxJacDPos,
				    flux_jac_type & fluxJacUNeg,
				    flux_jac_type & fluxJacUPos,
				    edge_rec_type & uMinusHalfNeg,
				    edge_rec_type & uMinusHalfPos,
				    edge_rec_type & uPlusHalfNeg,
				    edge_rec_type & uPlusHalfPos) const
  {
    // for inner cells, velocity and Jacobian
    // are both computed according to the order selected by the user
    // because for inner cells we support also Jacobians for Weno

    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    constexpr int zAxis = 3;

    // allocate gradients of reconstructed states
    // the size depends on the scheme selected
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    reconstruction_gradient_t gradLNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradLPos(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRPos(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradBNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradBPos(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradFNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradFPos(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradDNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradDPos(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradUNeg(numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradUPos(numDofPerCell, stencilSize-1);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impl::ComputeDirectionalFluxBalanceJacobianOnInteriorCell<
	  pda::ee::impl::ComputeDirectionalFluxValuesAndJacobians<
	    pda::impl::ReconstructorForDiscreteFunction<
	      dimensionality, numDofPerCell, MeshType, U_t, edge_rec_type, reconstruction_gradient_t>,
	    numDofPerCell, scalar_type, flux_type, flux_jac_type>,
	  dimensionality, numDofPerCell, MeshType, jacobian_type>,
      numDofPerCell, V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.dxInv(),
		    /* end args for velo */
		    J, xAxis, m_meshObj,
		    /* end args for jac */
		    m_fluxEn, normalX_, m_gamma, fluxL, fluxR,
		    fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
		    /* end args for flux */
		    xAxis, toReconstructionScheme(m_recEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		    gradLNeg, gradLPos, gradRNeg, gradRPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.dyInv(),
		    /* end args for velo */
		    J, yAxis, m_meshObj,
		    /* end args for jac */
		    m_fluxEn, normalY_, m_gamma, fluxB, fluxF,
		    fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
		    /* end args for flux */
		    yAxis, toReconstructionScheme(m_recEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		    gradBNeg, gradBPos, gradFNeg, gradFPos
		    /* end args for reconstructor */
		    );

    functor_type Fz(V, m_meshObj.dzInv(),
		    /* end args for velo */
		    J, zAxis, m_meshObj,
		    /* end args for jac */
		    m_fluxEn, normalZ_, m_gamma, fluxD, fluxU,
		    fluxJacDNeg, fluxJacDPos, fluxJacUNeg, fluxJacUPos,
		    /* end args for flux */
		    zAxis, toReconstructionScheme(m_recEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		    gradDNeg, gradDPos, gradUNeg, gradUPos
		    /* end args for reconstructor */
		    );

    const auto & graphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      Fx(smPt);
      Fy(smPt);
      Fz(smPt);
    }
  }

  template<class U_t, class V_t>
  void velocityAndJacNearBDCellsImplFirstOrder(const U_t & U,
					       const scalar_type currentTime,
					       V_t & V,
					       jacobian_type & J,
					       flux_type & fluxL,
					       flux_type & fluxR,
					       flux_type & fluxB,
					       flux_type & fluxF,
					       flux_type & fluxD,
					       flux_type & fluxU,
					       flux_jac_type & fluxJacLNeg,
					       flux_jac_type & fluxJacLPos,
					       flux_jac_type & fluxJacRNeg,
					       flux_jac_type & fluxJacRPos,
					       flux_jac_type & fluxJacBNeg,
					       flux_jac_type & fluxJacBPos,
					       flux_jac_type & fluxJacFNeg,
					       flux_jac_type & fluxJacFPos,
					       flux_jac_type & fluxJacDNeg,
					       flux_jac_type & fluxJacDPos,
					       flux_jac_type & fluxJacUNeg,
					       flux_jac_type & fluxJacUPos,
					       edge_rec_type & uMinusHalfNeg,
					       edge_rec_type & uMinusHalfPos,
					       edge_rec_type & uPlusHalfNeg,
					       edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    constexpr int zAxis = 3;
    assert(m_recEn == InviscidFluxReconstruction::FirstOrder);

    // if here, then the scheme for velocity matches
    // the one for Jacobian so we can use same functors

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_container_type,
      U_t, MeshType, ghost_container_type>;
    stencil_filler_t FillStencilX(reconstructionTypeToStencilSize(m_recEn),
				   U, m_meshObj, m_ghostLeft, m_ghostRight,
				   m_stencilVals, xAxis);

    stencil_filler_t FillStencilY(reconstructionTypeToStencilSize(m_recEn),
				   U, m_meshObj, m_ghostBack, m_ghostFront,
				   m_stencilVals, yAxis);

    stencil_filler_t FillStencilZ(reconstructionTypeToStencilSize(m_recEn),
				   U, m_meshObj, m_ghostDown, m_ghostUp,
				   m_stencilVals, zAxis);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impl::ComputeDirectionalFluxBalanceFirstOrderJacobianOnBoundaryCell<
	  pda::ee::impl::ComputeDirectionalFluxValuesAndJacobians<
	    pda::impl::ReconstructorFromStencil<
	      numDofPerCell, edge_rec_type, stencil_container_type>,
	    numDofPerCell, scalar_type, flux_type, flux_jac_type>,
	  dimensionality, numDofPerCell, MeshType, jacobian_type>,
      numDofPerCell, V_t, scalar_type
      >;

    functor_type funcx(V, m_meshObj.dxInv(),
		       /* end args for velo */
		       J, xAxis, m_meshObj,
		       /* end args for jac */
		       m_fluxEn, normalX_, m_gamma, fluxL, fluxR,
		       fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
		       /* end args for flux */
		       toReconstructionScheme(m_recEn), m_stencilVals,
		       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		       /* end args for reconstructor */
		       );

    functor_type funcy(V, m_meshObj.dyInv(),
		       /* end args for velo */
		       J, yAxis, m_meshObj,
		       /* end args for jac */
		       m_fluxEn, normalY_, m_gamma, fluxB, fluxF,
		       fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
		       /* end args for flux */
		       toReconstructionScheme(m_recEn), m_stencilVals,
		       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		       /* end args for reconstructor */
		       );

    functor_type funcz(V, m_meshObj.dzInv(),
		       /* end args for velo */
		       J, zAxis, m_meshObj,
		       /* end args for jac */
		       m_fluxEn, normalZ_, m_gamma, fluxD, fluxU,
		       fluxJacDNeg, fluxJacDPos, fluxJacUNeg, fluxJacUPos,
		       /* end args for flux */
		       toReconstructionScheme(m_recEn), m_stencilVals,
		       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		       /* end args for reconstructor */
		       );

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

    const auto & graphRows = m_meshObj.graphRowsOfCellsNearBd();
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];

      FillStencilX(smPt, it);
      auto bcTypeX = findCellBdType(smPt, xAxis);
      const auto & factorsX = (bcTypeX == 1) ? bcCellJacFactorsReflectiveX : bcCellJacFactorsDefault;
      funcx(smPt, factorsX, bcTypeX);

      FillStencilY(smPt, it);
      auto bcTypeY = findCellBdType(smPt, yAxis);
      const auto & factorsY = (bcTypeY == 1) ? bcCellJacFactorsReflectiveY : bcCellJacFactorsDefault;
      funcy(smPt, factorsY, bcTypeY);

      FillStencilZ(smPt, it);
      auto bcTypeZ = findCellBdType(smPt, zAxis);
      const auto & factorsZ = (bcTypeZ == 1) ? bcCellJacFactorsReflectiveZ : bcCellJacFactorsDefault;
      funcz(smPt, factorsZ, bcTypeZ);
    }
  }

  template<class U_t, class V_t>
  void velocityAndJacNearBDCellsImplDifferentScheme(const U_t & U,
						    const scalar_type currentTime,
						    V_t & V,
						    jacobian_type & J,
						    flux_type & fluxL,
						    flux_type & fluxR,
						    flux_type & fluxB,
						    flux_type & fluxF,
						    flux_type & fluxD,
						    flux_type & fluxU,
						    flux_jac_type & fluxJacLNeg,
						    flux_jac_type & fluxJacLPos,
						    flux_jac_type & fluxJacRNeg,
						    flux_jac_type & fluxJacRPos,
						    flux_jac_type & fluxJacBNeg,
						    flux_jac_type & fluxJacBPos,
						    flux_jac_type & fluxJacFNeg,
						    flux_jac_type & fluxJacFPos,
						    flux_jac_type & fluxJacDNeg,
						    flux_jac_type & fluxJacDPos,
						    flux_jac_type & fluxJacUNeg,
						    flux_jac_type & fluxJacUPos,
						    edge_rec_type & uMinusHalfNeg,
						    edge_rec_type & uMinusHalfPos,
						    edge_rec_type & uPlusHalfNeg,
						    edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    constexpr int zAxis = 3;

    // if here, then the velocity must be computed with Weno,
    /// while the jacobian must be computed with first order

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_container_type,
      U_t, MeshType, ghost_container_type>;

    // *****************************
    // *** functors for velocity ***
    // *****************************
    stencil_filler_t FillStencilVeloX(reconstructionTypeToStencilSize(m_recEn),
				      U, m_meshObj, m_ghostLeft, m_ghostRight,
				      m_stencilVals, xAxis);
    stencil_filler_t FillStencilVeloY(reconstructionTypeToStencilSize(m_recEn),
				      U, m_meshObj, m_ghostBack, m_ghostFront,
				      m_stencilVals, yAxis);
    stencil_filler_t FillStencilVeloZ(reconstructionTypeToStencilSize(m_recEn),
				      U, m_meshObj, m_ghostDown, m_ghostUp,
				      m_stencilVals, zAxis);

    using velo_functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::ee::impl::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorFromStencil<
	    numDofPerCell, edge_rec_type, stencil_container_type>,
	  numDofPerCell, scalar_type, flux_type>,
      numDofPerCell, V_t, scalar_type
      >;

    velo_functor_type funcVeloX(V, m_meshObj.dxInv(),
				/* end args for velo */
				m_fluxEn, normalX_, m_gamma, fluxL, fluxR,
				/* end args for flux */
				toReconstructionScheme(m_recEn), m_stencilVals,
				uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
				/* end args for reconstructor */
				);

    velo_functor_type funcVeloY(V, m_meshObj.dyInv(),
				/* end args for velo */
				m_fluxEn, normalY_, m_gamma, fluxB, fluxF,
				/* end args for flux */
				toReconstructionScheme(m_recEn), m_stencilVals,
				uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
				/* end args for reconstructor */
				);

    velo_functor_type funcVeloZ(V, m_meshObj.dzInv(),
				/* end args for velo */
				m_fluxEn, normalZ_, m_gamma, fluxD, fluxU,
				/* end args for flux */
				toReconstructionScheme(m_recEn), m_stencilVals,
				uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
				/* end args for reconstructor */
				);

    // *****************************
    // *** functors for jacobian ***
    // *****************************
    const auto firstOrderRec = pda::InviscidFluxReconstruction::FirstOrder;
    const auto stencilSizeForJ = reconstructionTypeToStencilSize(firstOrderRec);
    stencil_container_type stencilValsForJ(numDofPerCell*stencilSizeForJ);
    stencil_filler_t FillStencilJacX(stencilSizeForJ,
				     U, m_meshObj, m_ghostLeft, m_ghostRight,
				     stencilValsForJ, xAxis);
    stencil_filler_t FillStencilJacY(stencilSizeForJ,
				     U, m_meshObj, m_ghostBack, m_ghostFront,
				     stencilValsForJ, yAxis);
    stencil_filler_t FillStencilJacZ(stencilSizeForJ,
				     U, m_meshObj, m_ghostDown, m_ghostUp,
				     stencilValsForJ, zAxis);

    using jac_functor_type =
      pda::impl::ComputeDirectionalFluxBalanceFirstOrderJacobianOnBoundaryCell<
	pda::ee::impl::ComputeDirectionalFluxJacobians<
	  pda::impl::ReconstructorFromStencil<
	    numDofPerCell, edge_rec_type, stencil_container_type>,
	  numDofPerCell, scalar_type, flux_jac_type>,
      dimensionality, numDofPerCell, MeshType, jacobian_type
      >;

    jac_functor_type funcJacX(J, xAxis, m_meshObj,
			      /* end args for jac */
			      m_fluxEn, normalX_, m_gamma,
			      fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
			      /* end args for flux */
			      toReconstructionScheme(firstOrderRec), stencilValsForJ,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			      /* end args for reconstructor */
			      );

    jac_functor_type funcJacY(J, yAxis, m_meshObj,
			      /* end args for jac */
			      m_fluxEn, normalY_, m_gamma,
			      fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
			      /* end args for flux */
			      toReconstructionScheme(firstOrderRec), stencilValsForJ,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			      /* end args for reconstructor */
			      );

    jac_functor_type funcJacZ(J, zAxis, m_meshObj,
			      /* end args for jac */
			      m_fluxEn, normalZ_, m_gamma,
			      fluxJacDNeg, fluxJacDPos, fluxJacUNeg, fluxJacUPos,
			      /* end args for flux */
			      toReconstructionScheme(firstOrderRec), stencilValsForJ,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			      /* end args for reconstructor */
			      );

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

    // ************
    // loop
    // ************
    const auto & graphRows = m_meshObj.graphRowsOfCellsNearBd();
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      FillStencilVeloX(smPt, it);
      funcVeloX(smPt);
      FillStencilJacX(smPt, it);
      auto bcTypeX = findCellBdType(smPt, xAxis);
      const auto & factorsX = (bcTypeX == 1) ? bcCellJacFactorsReflectiveX : bcCellJacFactorsDefault;
      funcJacX(smPt, factorsX, 1);

      FillStencilVeloY(smPt, it);
      funcVeloY(smPt);
      FillStencilJacY(smPt, it);
      auto bcTypeY = findCellBdType(smPt, yAxis);
      const auto & factorsY = (bcTypeY == 1) ? bcCellJacFactorsReflectiveY : bcCellJacFactorsDefault;
      funcJacY(smPt, factorsY, 1);

      FillStencilVeloZ(smPt, it);
      funcVeloZ(smPt);
      FillStencilJacZ(smPt, it);
      auto bcTypeZ = findCellBdType(smPt, zAxis);
      const auto & factorsZ = (bcTypeZ == 1) ? bcCellJacFactorsReflectiveZ : bcCellJacFactorsDefault;
      funcJacZ(smPt, factorsZ, 1);

    }
  }

  template<class U_t, class V_t>
  void velocityOnlyInnerCellsImpl(const U_t & U,
				  const scalar_type currentTime,
				  V_t & V,
				  flux_type & fluxL,
				  flux_type & fluxR,
				  flux_type & fluxB,
				  flux_type & fluxF,
				  flux_type & fluxD,
				  flux_type & fluxU,
				  edge_rec_type & uMinusHalfNeg,
				  edge_rec_type & uMinusHalfPos,
				  edge_rec_type & uPlusHalfNeg,
				  edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    constexpr int zAxis = 3;

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::ee::impl::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorForDiscreteFunction<
	    dimensionality, numDofPerCell, MeshType, U_t, edge_rec_type>,
	  numDofPerCell, scalar_type, flux_type>,
      numDofPerCell, V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.dxInv(),
		    /* end args for velo */
		    m_fluxEn, normalX_, m_gamma, fluxL, fluxR,
		    /* end args for flux */
		    xAxis, toReconstructionScheme(m_recEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.dyInv(),
		    /* end args for velo */
		    m_fluxEn, normalY_, m_gamma, fluxB, fluxF,
		    /* end args for flux */
		    yAxis, toReconstructionScheme(m_recEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fz(V, m_meshObj.dzInv(),
		    /* end args for velo */
		    m_fluxEn, normalZ_, m_gamma, fluxB, fluxF,
		    /* end args for flux */
		    zAxis, toReconstructionScheme(m_recEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		    /* end args for reconstructor */
		    );

    const auto & graphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      Fx(smPt);
      Fy(smPt);
      Fz(smPt);
    }
  }

  template<class U_t, class V_t>
  void velocityOnlyNearBdCellsImpl(const U_t & U,
				   const scalar_type currentTime,
				   V_t & V,
				   flux_type & fluxL,
				   flux_type & fluxR,
				   flux_type & fluxB,
				   flux_type & fluxF,
				   flux_type & fluxD,
				   flux_type & fluxU,
				   edge_rec_type & uMinusHalfNeg,
				   edge_rec_type & uMinusHalfPos,
				   edge_rec_type & uPlusHalfNeg,
				   edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    constexpr int zAxis = 3;

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_container_type,
      U_t, MeshType, ghost_container_type>;
    stencil_filler_t FillStencilX(reconstructionTypeToStencilSize(m_recEn),
				  U, m_meshObj, m_ghostLeft, m_ghostRight,
				  m_stencilVals, xAxis);

    stencil_filler_t FillStencilY(reconstructionTypeToStencilSize(m_recEn),
				  U, m_meshObj, m_ghostBack, m_ghostFront,
				  m_stencilVals, yAxis);

    stencil_filler_t FillStencilZ(reconstructionTypeToStencilSize(m_recEn),
				  U, m_meshObj, m_ghostDown, m_ghostUp,
				  m_stencilVals, zAxis);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::ee::impl::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorFromStencil<
	    numDofPerCell, edge_rec_type, stencil_container_type>,
	  numDofPerCell, scalar_type, flux_type>,
      numDofPerCell, V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.dxInv(),
		    /* end args for velo */
		    m_fluxEn, normalX_, m_gamma, fluxL, fluxR,
		    /* end args for flux */
		    toReconstructionScheme(m_recEn), m_stencilVals,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.dyInv(),
		    /* end args for velo */
		    m_fluxEn, normalY_, m_gamma, fluxB, fluxF,
		    /* end args for flux */
		    toReconstructionScheme(m_recEn), m_stencilVals,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fz(V, m_meshObj.dzInv(),
		    /* end args for velo */
		    m_fluxEn, normalZ_, m_gamma, fluxD, fluxU,
		    /* end args for flux */
		    toReconstructionScheme(m_recEn), m_stencilVals,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		    /* end args for reconstructor */
		    );

    const auto & graphRows = m_meshObj.graphRowsOfCellsNearBd();
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      FillStencilX(smPt, it);
      Fx(smPt);
      FillStencilY(smPt, it);
      Fy(smPt);
      FillStencilZ(smPt, it);
      Fz(smPt);
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
    const auto stencilSize    = reconstructionTypeToStencilSize(m_recEn);
    const auto numGhostValues = numDofPerCell*((stencilSize-1)/2);

    const index_t s1 = m_meshObj.numCellsBd();
    ::pressiodemoapps::resize(m_ghostLeft,  s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostFront, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostBack,  s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostDown, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostUp,  s1, numGhostValues);
  }

protected:
  scalar_type m_gamma = static_cast<scalar_type>(1.4);

  const MeshType & m_meshObj;
  ::pressiodemoapps::Euler3d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  ::pressiodemoapps::InviscidFluxScheme m_fluxEn;

  // note that dof refers to the degress of freedom,
  // which is NOT same as grid points.
  // SampleMesh_ identifies the velocity/residual locations
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  mutable stencil_container_type m_stencilVals;

  mutable ghost_container_type m_ghostLeft;
  mutable ghost_container_type m_ghostRight;
  mutable ghost_container_type m_ghostBack;
  mutable ghost_container_type m_ghostFront;
  mutable ghost_container_type m_ghostDown;
  mutable ghost_container_type m_ghostUp;

  const std::array<scalar_type, dimensionality> normalX_{1, 0, 0};
  const std::array<scalar_type, dimensionality> normalY_{0, 1, 0};
  const std::array<scalar_type, dimensionality> normalZ_{0, 0, 1};
};

}}}//end namespace pressiodemoapps::ee::impl
#endif









//   scalar_type gamma()           const{ return m_gamma; }
//   index_t totalDofSampleMesh()  const{ return m_numDofSampleMesh; }
//   index_t totalDofStencilMesh() const{ return m_numDofStencilMesh; }

//   velocity_type createVelocity() const
//   {
//     velocity_type V(m_numDofSampleMesh);
//     return V;
//   }

//   void velocity(const state_type & state,
// 		const scalar_type timeValue,
// 		velocity_type & V) const
//   {
//     flux_t FL(numDofPerCell);
//     flux_t FR(numDofPerCell);
//     flux_t FB(numDofPerCell);
//     flux_t FF(numDofPerCell);
//     flux_t FD(numDofPerCell);
//     flux_t FU(numDofPerCell);
//     edge_rec_t uMinusHalfNeg(numDofPerCell);
//     edge_rec_t uMinusHalfPos(numDofPerCell);
//     edge_rec_t uPlusHalfNeg (numDofPerCell);
//     edge_rec_t uPlusHalfPos (numDofPerCell);

//     fillGhosts(state);

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//     if (!m_onlyComputeVelocity)
//     {
//       auto values = m_jacobian.valuePtr();
//       for (int i=0; i<m_jacobian.nonZeros(); ++i){
// 	values[i] = 0;
//       }
//     }
// #endif

//     velocityCellsNearBdImpl(state, timeValue,
// 			    V, FL, FR, FB, FF, FD, FU,
// 			    uMinusHalfNeg, uMinusHalfPos,
// 			    uPlusHalfNeg,  uPlusHalfPos);
//     velocityInnerCellsImpl(state, timeValue,
// 			   V, FL, FR, FB, FF, FD, FU,
// 			   uMinusHalfNeg, uMinusHalfPos,
// 			   uPlusHalfNeg,  uPlusHalfPos);
//   }

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//   jacobian_type createJacobian() const{
//     jacobian_type JJ(m_numDofSampleMesh, m_numDofStencilMesh);
//     return JJ;
//   }

//   void jacobian(const state_type & state,
// 		const scalar_type time,
// 		jacobian_type & J) const
//   {
//     if (!m_onlyComputeVelocity){
//       // relies on jacobian been computed in velocity
//       J = m_jacobian;
//     }
//   }

//   // the Jacobian is by default fused with the velocity,
//   // this method allows one to disable the jacobian
//   // so only velocity is computed
//   void disableJacobian() {
//     m_onlyComputeVelocity = true;
//   }
// #endif

// private:
//   void computeDofs(){
//     // calculate total num of dofs on sample and stencil mesh
//     m_numDofStencilMesh = m_meshObj.stencilMeshSize() * numDofPerCell;
//     m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * numDofPerCell;
//   }

//   void allocateStencilValuesContainer()
//   {
//     // the stencil size needed is determined by the desired reconstruction
//     // kind NOT from the mesh. THis is important because for example
//     // the mesh can have a wider connectivity that what is needed.
//     const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_recEn);
//     ::pressiodemoapps::resize(m_stencilVals, numDofPerCell*stencilSizeNeeded);
//   }

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//   void initializeJacobian()
//   {
//     initializeJacobianFirstOrder();

//     if (!m_jacobian.isCompressed()){
//       m_jacobian.makeCompressed();
//     }

//     m_jacNonZerosCount = m_jacobian.nonZeros();

//     // if Jacobian is disabled, free it
//     if (m_onlyComputeVelocity){
//       ::pressiodemoapps::resize(m_jacobian, 0, 0);
//     }
//   }

//   void initializeJacobianFirstOrder()
//   {
//     using Tr = Eigen::Triplet<scalar_type>;
//     std::vector<Tr> trList;

//     const scalar_type val0 = 0.;
//     const auto & graph = m_meshObj.graph();
//     for (int cell=0; cell<m_meshObj.sampleMeshSize(); ++cell)
//       {
// 	const auto jacRowOfCurrentCellDensity = cell*numDofPerCell;
// 	const auto ci0  = graph(cell, 0)*numDofPerCell;

// 	for (int k=0; k<numDofPerCell; ++k){
// 	  for (int j=0; j<numDofPerCell; ++j){
// 	    trList.push_back( Tr(jacRowOfCurrentCellDensity+k, ci0+j, val0) );
// 	  }
// 	}

// 	const auto L0 = graph(cell, 1);
// 	const auto F0 = graph(cell, 2);
// 	const auto R0 = graph(cell, 3);
// 	const auto B0 = graph(cell, 4);
// 	const auto D0 = graph(cell, 5);
// 	const auto U0 = graph(cell, 6);

// 	for (int i=1; i<=6; ++i)
// 	{
// 	  const auto gID = graph(cell, i);

// 	  if (gID != -1){
// 	    const auto ci = gID*numDofPerCell;
// 	    for (int k=0; k<numDofPerCell; ++k){
// 	      for (int j=0; j<numDofPerCell; ++j){
// 		trList.push_back( Tr(jacRowOfCurrentCellDensity+k, ci+j, val0) );
// 	      }
// 	    }
// 	  }
// 	}
//       }

//     m_jacobian.setFromTriplets(trList.begin(), trList.end());
//   }
// #endif

//   void velocityInnerCellsImpl(const state_type & U,
// 			      const scalar_type timeIn,
// 			      velocity_type & V,
// 			      flux_t & FL,
// 			      flux_t & FR,
// 			      flux_t & FB,
// 			      flux_t & FF,
// 			      flux_t & FD,
// 			      flux_t & FU,
// 			      edge_rec_t & uMinusHalfNeg,
// 			      edge_rec_t & uMinusHalfPos,
// 			      edge_rec_t & uPlusHalfNeg,
// 			      edge_rec_t & uPlusHalfPos) const
//   {

//     constexpr int xAxis = 1;
//     constexpr int yAxis = 2;
//     constexpr int zAxis = 3;
//     const auto dxInv   = m_meshObj.dxInv();
//     const auto dyInv   = m_meshObj.dyInv();
//     const auto dzInv   = m_meshObj.dzInv();

//     using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromState<
//       dimensionality, edge_rec_t, state_type, mesh_t>;

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//     using jac_fnct_t = ::pressiodemoapps::impl::FirstOrderInnerCellJacobianFunctor<
//       dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, mesh_t>;
// #endif

//     // reconstruct functor for face fluxes
//     // here we need to use whatever order (m_recEn) user decides
//     rec_fnct_t ReconstructorX(xAxis, m_recEn, U, m_meshObj,
// 			      uMinusHalfNeg, uMinusHalfPos,
// 			      uPlusHalfNeg,  uPlusHalfPos);

//     rec_fnct_t ReconstructorY(yAxis, m_recEn, U, m_meshObj,
// 			      uMinusHalfNeg, uMinusHalfPos,
// 			      uPlusHalfNeg,  uPlusHalfPos);

//     rec_fnct_t ReconstructorZ(zAxis, m_recEn, U, m_meshObj,
// 			      uMinusHalfNeg, uMinusHalfPos,
// 			      uPlusHalfNeg,  uPlusHalfPos);

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//     // cell jacobian functor
//     // currently, REGARDLESS of the reconstruction scheme,
//     // we have only first-order Jacobian so we need a first-order reconstructor for the Jacobian
//     flux_jac_type JLneg, JLpos;
//     flux_jac_type JRneg, JRpos;
//     flux_jac_type JBneg, JBpos;
//     flux_jac_type JFneg, JFpos;
//     flux_jac_type JDneg, JDpos;
//     flux_jac_type JUneg, JUpos;

//     const auto firstOrder = ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder;
//     edge_rec_t uMinusHalfNegForJ(numDofPerCell);
//     edge_rec_t uMinusHalfPosForJ(numDofPerCell);
//     edge_rec_t uPlusHalfNegForJ(numDofPerCell);
//     edge_rec_t uPlusHalfPosForJ(numDofPerCell);
//     rec_fnct_t ReconstructorXForJ(xAxis, firstOrder, U, m_meshObj,
// 				  uMinusHalfNegForJ, uMinusHalfPosForJ,
// 				  uPlusHalfNegForJ,  uPlusHalfPosForJ);

//     rec_fnct_t ReconstructorYForJ(yAxis, firstOrder, U, m_meshObj,
// 				  uMinusHalfNegForJ, uMinusHalfPosForJ,
// 				  uPlusHalfNegForJ,  uPlusHalfPosForJ);

//     rec_fnct_t ReconstructorZForJ(zAxis, firstOrder, U, m_meshObj,
// 				  uMinusHalfNegForJ, uMinusHalfPosForJ,
// 				  uPlusHalfNegForJ,  uPlusHalfPosForJ);

//     jac_fnct_t CellJacobianFunctorX(m_jacobian, m_meshObj, JLneg, JLpos, JRneg, JRpos, xAxis);
//     jac_fnct_t CellJacobianFunctorY(m_jacobian, m_meshObj, JBneg, JBpos, JFneg, JFpos, yAxis);
//     jac_fnct_t CellJacobianFunctorZ(m_jacobian, m_meshObj, JDneg, JDpos, JUneg, JUpos, zAxis);
// #endif

//     // deal with cells away from boundaries
//     const auto & rowsIn = m_meshObj.graphRowsOfCellsAwayFromBd();
//     for (int it=0; it<rowsIn.size(); ++it)
//     {
//       const auto smPt = rowsIn[it];
//       const auto vIndex = smPt*numDofPerCell;

//       // *** X ***
//       ReconstructorX.template operator()<numDofPerCell>(smPt);

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//       if (!m_onlyComputeVelocity){
// 	// note that REGARDLESS of the reconstruction scheme,
// 	// we currently only have only first-order Jacobian so we need
// 	// to run the reconstructor for the Jacobian
// 	// which will ensure that uMinusNegForJ, etc have the right values
// 	ReconstructorXForJ.template operator()<numDofPerCell>(smPt);
//       }
// #endif

//       switch(m_fluxEn)
// 	{
// 	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
// 	  ee_rusanov_flux_five_dof(FL, uMinusHalfNeg, uMinusHalfPos, normalX_, m_gamma);
// 	  ee_rusanov_flux_five_dof(FR, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gamma);

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
// 	  if (!m_onlyComputeVelocity){
// 	    ee_rusanov_flux_jacobian_five_dof(JLneg, JLpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
// 					 normalX_, m_gamma);
// 	    ee_rusanov_flux_jacobian_five_dof(JRneg, JRpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
// 					 normalX_, m_gamma);
// 	  }
// #endif
// 	  break;
// 	}

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//       if (!m_onlyComputeVelocity){
// 	CellJacobianFunctorX(smPt);
//       }
// #endif

//       // *** Y ***
//       ReconstructorY.template operator()<numDofPerCell>(smPt);

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//       if (!m_onlyComputeVelocity){
// 	// note that REGARDLESS of the reconstruction scheme,
// 	// we currently only have only first-order Jacobian so we need
// 	// to run the reconstructor for the Jacobian
// 	// which will ensure that uMinusNegForJ, etc have the right values
// 	ReconstructorYForJ.template operator()<numDofPerCell>(smPt);
//       }
// #endif

//       switch(m_fluxEn)
// 	{
// 	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
// 	  ee_rusanov_flux_five_dof(FB, uMinusHalfNeg, uMinusHalfPos, normalY_, m_gamma);
// 	  ee_rusanov_flux_five_dof(FF, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gamma);

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
// 	  if (!m_onlyComputeVelocity){
// 	    ee_rusanov_flux_jacobian_five_dof(JBneg, JBpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
// 					 normalY_, m_gamma);
// 	    ee_rusanov_flux_jacobian_five_dof(JFneg, JFpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
// 					 normalY_, m_gamma);
// 	  }
// #endif
// 	  break;
// 	}

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//       if (!m_onlyComputeVelocity){
// 	CellJacobianFunctorY(smPt);
//       }
// #endif

//       // *** Z ***
//       ReconstructorZ.template operator()<numDofPerCell>(smPt);

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//       if (!m_onlyComputeVelocity){
// 	// note that REGARDLESS of the reconstruction scheme,
// 	// we currently only have only first-order Jacobian so we need
// 	// to run the reconstructor for the Jacobian
// 	// which will ensure that uMinusNegForJ, etc have the right values
// 	ReconstructorZForJ.template operator()<numDofPerCell>(smPt);
//       }
// #endif

//       switch(m_fluxEn)
// 	{
// 	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
// 	  ee_rusanov_flux_five_dof(FD, uMinusHalfNeg, uMinusHalfPos, normalZ_, m_gamma);
// 	  ee_rusanov_flux_five_dof(FU, uPlusHalfNeg,  uPlusHalfPos,  normalZ_, m_gamma);

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
// 	  if (!m_onlyComputeVelocity){
// 	    ee_rusanov_flux_jacobian_five_dof(JDneg, JDpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
// 					 normalZ_, m_gamma);
// 	    ee_rusanov_flux_jacobian_five_dof(JUneg, JUpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
// 					 normalZ_, m_gamma);
// 	  }
// #endif
// 	  break;
// 	}

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//       if (!m_onlyComputeVelocity){
// 	CellJacobianFunctorZ(smPt);
//       }
// #endif

//       V(vIndex)   = dxInv*(FL(0) - FR(0)) + dyInv*(FB(0) - FF(0)) + dzInv*(FD(0) - FU(0));
//       V(vIndex+1) = dxInv*(FL(1) - FR(1)) + dyInv*(FB(1) - FF(1)) + dzInv*(FD(1) - FU(1));
//       V(vIndex+2) = dxInv*(FL(2) - FR(2)) + dyInv*(FB(2) - FF(2)) + dzInv*(FD(2) - FU(2));
//       V(vIndex+3) = dxInv*(FL(3) - FR(3)) + dyInv*(FB(3) - FF(3)) + dzInv*(FD(3) - FU(3));
//       V(vIndex+4) = dxInv*(FL(4) - FR(4)) + dyInv*(FB(4) - FF(4)) + dzInv*(FD(4) - FU(4));
//     }
//   }

//   int findCellBdType(index_t graphRow, int axis) const
//   {
//     // 0: Neumann
//     // 1: Reflective
//     // 2: Dirichlet
//     constexpr int neumann = 0;
//     constexpr int reflective = 1;
//     constexpr int dirichlet = 2;

//     if (m_probEn == ::pressiodemoapps::Euler3d::SedovSymmetry)
//     {
//       if (axis == 1 && m_meshObj.hasBdLeft3d(graphRow)){
// 	return reflective;
//       }

//       if (axis == 2 && m_meshObj.hasBdBack3d(graphRow)){
// 	return reflective;
//       }

//       if (axis == 3 && m_meshObj.hasBdBottom3d(graphRow)){
// 	return reflective;
//       }

//       return neumann;
//     }

//     return 0;
//   }

//   void velocityCellsNearBdImpl(const state_type & U,
// 			       const scalar_type timeIn,
// 			       velocity_type & V,
// 			       flux_t & FL,
// 			       flux_t & FR,
// 			       flux_t & FB,
// 			       flux_t & FF,
// 			       flux_t & FD,
// 			       flux_t & FU,
// 			       edge_rec_t & uMinusHalfNeg,
// 			       edge_rec_t & uMinusHalfPos,
// 			       edge_rec_t & uPlusHalfNeg,
// 			       edge_rec_t & uPlusHalfPos) const
//   {

//     constexpr int xAxis = 1;
//     constexpr int yAxis = 2;
//     constexpr int zAxis = 3;

//     const auto dxInv   = m_meshObj.dxInv();
//     const auto dyInv   = m_meshObj.dyInv();
//     const auto dzInv   = m_meshObj.dzInv();

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//     flux_jac_type JLneg, JLpos;
//     flux_jac_type JRneg, JRpos;
//     flux_jac_type JBneg, JBpos;
//     flux_jac_type JFneg, JFpos;
//     flux_jac_type JDneg, JDpos;
//     flux_jac_type JUneg, JUpos;
// #endif

//     using sfiller_t = ::pressiodemoapps::impl::StencilFiller<
//       dimensionality, numDofPerCell, stencil_values_t, state_type, mesh_t, ghost_t>;

//     using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromStencilFiveDofPerCell<
//       edge_rec_t, stencil_values_t>;

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//     using jac_fnct_t = ::pressiodemoapps::impl::FirstOrderBdCellJacobianFunctor<
//       dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, mesh_t>;
// #endif

//     const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
//     sfiller_t StencilFillerX(stencilSize, U, m_meshObj,
// 			     m_ghostLeft, m_ghostRight,
// 			     m_stencilVals, xAxis);

//     sfiller_t StencilFillerY(stencilSize, U, m_meshObj,
// 			     m_ghostBack, m_ghostFront,
// 			     m_stencilVals, yAxis);

//     sfiller_t StencilFillerZ(stencilSize, U, m_meshObj,
// 			     m_ghostBottom, m_ghostTop,
// 			     m_stencilVals, zAxis);

//     rec_fnct_t Reconstructor(m_recEn, m_stencilVals,
// 			     uMinusHalfNeg, uMinusHalfPos,
// 			     uPlusHalfNeg,  uPlusHalfPos);


// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//     // ----------------------------------------
//     // functors needed to compute cell jacobian
//     // we have only first-order Jacobian for now, so we need
//     // dedicate functors because we cannot use those above.
//     // Once we have cell Jacobians of various order, we can change this.
//     // ----------------------------------------
//     const auto cellJacOrdEn = ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder;

//     stencil_values_t stencilValsForJ = {};
//     const auto stencilSizeForJ = reconstructionTypeToStencilSize(cellJacOrdEn);
//     ::pressiodemoapps::resize(stencilValsForJ, numDofPerCell*stencilSizeForJ);
//     sfiller_t FillStencilValuesXFunctorForJ(stencilSizeForJ, U, m_meshObj,
// 					    m_ghostLeft, m_ghostRight, stencilValsForJ, xAxis);

//     sfiller_t FillStencilValuesYFunctorForJ(stencilSizeForJ, U, m_meshObj,
// 					    m_ghostBack, m_ghostFront, stencilValsForJ, yAxis);

//     sfiller_t FillStencilValuesZFunctorForJ(stencilSizeForJ, U, m_meshObj,
// 					    m_ghostBottom, m_ghostTop, stencilValsForJ, zAxis);

//     edge_rec_t uMinusHalfNegForJ(numDofPerCell);
//     edge_rec_t uMinusHalfPosForJ(numDofPerCell);
//     edge_rec_t uPlusHalfNegForJ(numDofPerCell);
//     edge_rec_t uPlusHalfPosForJ(numDofPerCell);
//     rec_fnct_t FaceValuesReconstructFunctorForJ(cellJacOrdEn, stencilValsForJ,
// 						uMinusHalfNegForJ, uMinusHalfPosForJ,
// 						uPlusHalfNegForJ,  uPlusHalfPosForJ);

//     jac_fnct_t CellJacobianFunctorX(m_jacobian, m_meshObj, JLneg, JLpos, JRneg, JRpos, xAxis);
//     jac_fnct_t CellJacobianFunctorY(m_jacobian, m_meshObj, JBneg, JBpos, JFneg, JFpos, yAxis);
//     jac_fnct_t CellJacobianFunctorZ(m_jacobian, m_meshObj, JDneg, JDpos, JUneg, JUpos, zAxis);

//     std::array<scalar_type, numDofPerCell> bcCellJacFactorsDefault;
//     bcCellJacFactorsDefault.fill(static_cast<scalar_type>(1));

//     std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveX;
//     std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveY;
//     std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveZ;

//     bcCellJacFactorsReflectiveX.fill(static_cast<scalar_type>(1));
//     bcCellJacFactorsReflectiveY.fill(static_cast<scalar_type>(1));
//     bcCellJacFactorsReflectiveZ.fill(static_cast<scalar_type>(1));
//     bcCellJacFactorsReflectiveX[1] = static_cast<scalar_type>(-1);
//     bcCellJacFactorsReflectiveY[2] = static_cast<scalar_type>(-1);
//     bcCellJacFactorsReflectiveZ[3] = static_cast<scalar_type>(-1);
// #endif

//     // -----
//     // loop
//     // -----
//     const auto & specialRows = m_meshObj.graphRowsOfCellsNearBd();
//     for (std::size_t it=0; it<specialRows.size(); ++it)
//     {
//       const auto smPt = specialRows[it];

//       // ------------
//       // X
//       // ------------
//       StencilFillerX(smPt, it);
//       Reconstructor();

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//       if (!m_onlyComputeVelocity){
// 	FillStencilValuesXFunctorForJ(smPt, it);
// 	FaceValuesReconstructFunctorForJ();
//       }
// #endif

//       switch(m_fluxEn)
// 	{
// 	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
// 	  ee_rusanov_flux_five_dof(FL, uMinusHalfNeg, uMinusHalfPos, normalX_, m_gamma);
// 	  ee_rusanov_flux_five_dof(FR, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gamma);

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
// 	  if (!m_onlyComputeVelocity){
// 	    ee_rusanov_flux_jacobian_five_dof(JLneg, JLpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
// 					 normalX_, m_gamma);
// 	    ee_rusanov_flux_jacobian_five_dof(JRneg, JRpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
// 					 normalX_, m_gamma);
// 	  }
// #endif

// 	  break;
// 	}

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//       if (!m_onlyComputeVelocity){
// 	auto bcType = findCellBdType(smPt, xAxis);
// 	const auto & factorsX = (bcType == 1) ? bcCellJacFactorsReflectiveX : bcCellJacFactorsDefault;
// 	CellJacobianFunctorX(smPt, factorsX, bcType);
//       }
// #endif

//       // ------------
//       // Y
//       // ------------
//       StencilFillerY(smPt, it);
//       Reconstructor();

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//       if (!m_onlyComputeVelocity){
//       FillStencilValuesYFunctorForJ(smPt, it);
//       FaceValuesReconstructFunctorForJ();
//       }
// #endif

//       switch(m_fluxEn)
// 	{
// 	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
// 	  ee_rusanov_flux_five_dof(FB, uMinusHalfNeg, uMinusHalfPos, normalY_, m_gamma);
// 	  ee_rusanov_flux_five_dof(FF, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gamma);

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
// 	  if (!m_onlyComputeVelocity){
// 	  ee_rusanov_flux_jacobian_five_dof(JBneg, JBpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
// 				       normalY_, m_gamma);
// 	  ee_rusanov_flux_jacobian_five_dof(JFneg, JFpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
// 				       normalY_, m_gamma);
// 	  }
// #endif

// 	  break;
// 	}

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//       if (!m_onlyComputeVelocity){
// 	auto bcType = findCellBdType(smPt, yAxis);
// 	const auto & factorsY = (bcType == 1) ? bcCellJacFactorsReflectiveY : bcCellJacFactorsDefault;
// 	CellJacobianFunctorY(smPt, factorsY, bcType);
//       }
// #endif

//       // ------------
//       // Z
//       // ------------
//       StencilFillerZ(smPt, it);
//       Reconstructor();

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//       if (!m_onlyComputeVelocity){
// 	FillStencilValuesZFunctorForJ(smPt, it);
// 	FaceValuesReconstructFunctorForJ();
//       }
// #endif

//       switch(m_fluxEn)
// 	{
// 	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
// 	  ee_rusanov_flux_five_dof(FD, uMinusHalfNeg, uMinusHalfPos, normalZ_, m_gamma);
// 	  ee_rusanov_flux_five_dof(FU, uPlusHalfNeg,  uPlusHalfPos,  normalZ_, m_gamma);

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
// 	  if (!m_onlyComputeVelocity){
// 	  ee_rusanov_flux_jacobian_five_dof(JDneg, JDpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
// 				       normalZ_, m_gamma);
// 	  ee_rusanov_flux_jacobian_five_dof(JUneg, JUpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
// 				       normalZ_, m_gamma);
// 	  }
// #endif

// 	  break;
// 	}

// #ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
//       if (!m_onlyComputeVelocity){
// 	auto bcType = findCellBdType(smPt, zAxis);
// 	const auto & factorsZ = (bcType == 1) ? bcCellJacFactorsReflectiveZ : bcCellJacFactorsDefault;
// 	CellJacobianFunctorZ(smPt, factorsZ, bcType);
//       }
// #endif

//       const auto vIndex = smPt*numDofPerCell;
//       V(vIndex)   = dxInv*(FL(0) - FR(0)) + dyInv*(FB(0) - FF(0)) + dzInv*(FD(0) - FU(0));
//       V(vIndex+1) = dxInv*(FL(1) - FR(1)) + dyInv*(FB(1) - FF(1)) + dzInv*(FD(1) - FU(1));
//       V(vIndex+2) = dxInv*(FL(2) - FR(2)) + dyInv*(FB(2) - FF(2)) + dzInv*(FD(2) - FU(2));
//       V(vIndex+3) = dxInv*(FL(3) - FR(3)) + dyInv*(FB(3) - FF(3)) + dzInv*(FD(3) - FU(3));
//       V(vIndex+4) = dxInv*(FL(4) - FR(4)) + dyInv*(FB(4) - FF(4)) + dzInv*(FD(4) - FU(4));
//     }
//   }


//   void fillGhosts(const state_type & U) const
//   {
//     const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
//     if (m_probEn == ::pressiodemoapps::Euler3d::SedovSymmetry)
//     {
//       using ghost_filler_t  = ::pressiodemoapps::ee::impl::Ghost3dSedov<state_type, mesh_t, ghost_t>;
//       ghost_filler_t ghF(stencilSize,   U, m_meshObj,
// 			 m_ghostLeft,   m_ghostRight,
// 			 m_ghostBack,   m_ghostFront,
// 			 m_ghostBottom, m_ghostTop);

//       const auto & rowsBd = m_meshObj.graphRowsOfCellsNearBd();

//       if (stencilSize==3){
// 	for (int it=0; it<rowsBd.size(); ++it){
// 	  ghF.template operator()<3>(rowsBd[it], it);
// 	}
//       }
//       else if (stencilSize==5){
// 	for (int it=0; it<rowsBd.size(); ++it){
// 	  ghF.template operator()<5>(rowsBd[it], it);
// 	}
//       }else{
// 	throw std::runtime_error("missing impl");
//       }
//     }
//   }

// private:
//   void allocateGhosts()
//   {
//     /*
//       for stencil = 3, at leftBoundary:
//       ---------------
//       |	 0,1,2,3,4 ||
//       | rho,       ||
//       | rho u,	   ||
//       | rho v,     ||
//       | rho w,     ||
//       | E	   ||
//       ---------------

//       for stencil = 7, at leftBoundary:
//       --------------------------------------
//       |	10,11,12,13,14  | 5,6,7,8,9 |  0,1,2,3,4 ||
//       |	     	        |
//       | rho,            | rho,       | rho        ||
//       | rho u,	        | rho*u      | rho*u      ||
//       | rho v,          | rho*v      | rho*v      ||
//       | rho w,          | rho*w      | rho*w      ||
//       | E	        | E          | E          ||
//       ---------------------------------------
//      */

//     const auto stencilSize    = reconstructionTypeToStencilSize(m_recEn);
//     const auto numGhostValues = numDofPerCell*((stencilSize-1)/2);

//     const index_t s1 = m_meshObj.numCellsBd();
//     ::pressiodemoapps::resize(m_ghostLeft,   s1, numGhostValues);
//     ::pressiodemoapps::resize(m_ghostRight,  s1, numGhostValues);
//     ::pressiodemoapps::resize(m_ghostBack,   s1, numGhostValues);
//     ::pressiodemoapps::resize(m_ghostFront,  s1, numGhostValues);
//     ::pressiodemoapps::resize(m_ghostBottom, s1, numGhostValues);
//     ::pressiodemoapps::resize(m_ghostTop,    s1, numGhostValues);
//   }
