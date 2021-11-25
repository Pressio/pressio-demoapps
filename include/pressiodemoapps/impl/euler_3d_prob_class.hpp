
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

  static constexpr int dimensionality{3};
  static constexpr int numDofPerCell{5};

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
      funcJacX(smPt, factorsX, bcTypeX);

      FillStencilVeloY(smPt, it);
      funcVeloY(smPt);
      FillStencilJacY(smPt, it);
      auto bcTypeY = findCellBdType(smPt, yAxis);
      const auto & factorsY = (bcTypeY == 1) ? bcCellJacFactorsReflectiveY : bcCellJacFactorsDefault;
      funcJacY(smPt, factorsY, bcTypeY);

      FillStencilVeloZ(smPt, it);
      funcVeloZ(smPt);
      FillStencilJacZ(smPt, it);
      auto bcTypeZ = findCellBdType(smPt, zAxis);
      const auto & factorsZ = (bcTypeZ == 1) ? bcCellJacFactorsReflectiveZ : bcCellJacFactorsDefault;
      funcJacZ(smPt, factorsZ, bcTypeZ);

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

template<class MeshType> constexpr int EigenEuler3dApp<MeshType>::numDofPerCell;
template<class MeshType> constexpr int EigenEuler3dApp<MeshType>::dimensionality;

}}}//end namespace pressiodemoapps::ee::impl
#endif
