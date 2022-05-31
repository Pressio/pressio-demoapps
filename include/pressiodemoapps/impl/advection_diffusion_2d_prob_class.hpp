
#ifndef PRESSIODEMOAPPS_ADVECTION_DIFFUSION_2D_PROB_CLASS_HPP_
#define PRESSIODEMOAPPS_ADVECTION_DIFFUSION_2D_PROB_CLASS_HPP_

#include "advection_diffusion_burgers_ghost_filler.hpp"
#include "advection_diffusion_2d_flux_functions.hpp"
#include "advection_diffusion_2d_initial_condition.hpp"
#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "advection_diffusion_2d_flux_mixin.hpp"
#include "mixin_directional_flux_balance.hpp"
#include "mixin_directional_flux_balance_jacobian.hpp"
#include "Eigen/Sparse"

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#include <omp.h>
#endif

namespace pressiodemoapps{ namespace impladvdiff2d{

// tags are used inside he public create function: create_problem_...()
// in the file ../advection_diffusion.hpp
// to dispatch to the proper problem
// so add new ones if a new problem is added
struct TagBurgersPeriodic{};
struct TagBurgersDirichlet{};


/////////////////////////
// eigen class
/////////////////////////
template<class MeshType>
class EigenApp
{

public:
  // required public aliases
  using index_t	      = typename MeshType::index_t;
  using scalar_type   = typename MeshType::scalar_t;
  using state_type    = Eigen::Matrix<scalar_type,Eigen::Dynamic,1>;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, index_t>;

private:
  static constexpr int dimensionality{2};

  using ghost_container_type      = Eigen::Matrix<scalar_type,
						  Eigen::Dynamic,
						  Eigen::Dynamic,
						  Eigen::RowMajor>;
  using stencil_container_type    = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using flux_type	          = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using edge_rec_type	          = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using flux_jac_type             = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;
  using reconstruction_gradient_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;

public:

  //
  // constructor for Burgers2d periodic
  //
  EigenApp(TagBurgersPeriodic /*tag*/,
	   const MeshType & meshObj,
	   ::pressiodemoapps::InviscidFluxReconstruction inviscidFluxRecEn,
	   ::pressiodemoapps::InviscidFluxScheme invFluxSchemeEn,
	   ::pressiodemoapps::ViscousFluxReconstruction visFluxRecEn,
	   scalar_type icPulseMagnitude,
	   scalar_type icSpread,
	   scalar_type diffusionCoeff,
	   scalar_type x0,
	   scalar_type y0)
    : m_numDofPerCell(2),
      m_probEn(::pressiodemoapps::AdvectionDiffusion2d::BurgersPeriodic),
      m_inviscidFluxRecEn(inviscidFluxRecEn),
      m_inviscidFluxSchemeEn(invFluxSchemeEn),
      m_viscousFluxRecEn(visFluxRecEn),
      m_burgers2d_icPulse(icPulseMagnitude),
      m_burgers2d_icSpread(icSpread),
      m_burgers2d_diffusion(diffusionCoeff),
      m_burgers2d_x0(x0),
      m_burgers2d_y0(y0),
      m_meshObj(meshObj)
  {
    m_numDofStencilMesh = m_meshObj.stencilMeshSize() * m_numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * m_numDofPerCell;
    // don't need to allocate ghosts because it is periodic
  }

  //
  // constructor for Burgers2d homogen dirichlet
  //
  EigenApp(TagBurgersDirichlet /*tag*/,
	   const MeshType & meshObj,
	   ::pressiodemoapps::InviscidFluxReconstruction inviscidFluxRecEn,
	   ::pressiodemoapps::InviscidFluxScheme invFluxSchemeEn,
	   ::pressiodemoapps::ViscousFluxReconstruction visFluxRecEn,
	   scalar_type icPulseMagnitude,
	   scalar_type icSpread,
	   scalar_type diffusionCoeff,
	   scalar_type x0,
	   scalar_type y0)
    : m_numDofPerCell(2),
      m_probEn(::pressiodemoapps::AdvectionDiffusion2d::BurgersDirichlet),
      m_inviscidFluxRecEn(inviscidFluxRecEn),
      m_inviscidFluxSchemeEn(invFluxSchemeEn),
      m_viscousFluxRecEn(visFluxRecEn),
      m_burgers2d_icPulse(icPulseMagnitude),
      m_burgers2d_icSpread(icSpread),
      m_burgers2d_diffusion(diffusionCoeff),
      m_burgers2d_x0(x0),
      m_burgers2d_y0(y0),
      m_meshObj(meshObj)
  {
    m_numDofStencilMesh = m_meshObj.stencilMeshSize() * m_numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * m_numDofPerCell;
    allocateGhosts();
  }

  state_type initialCondition() const
  {
    state_type initialState(m_numDofStencilMesh);

    if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersDirichlet ||
	m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersPeriodic)
      {
	burgers2d_gaussian(initialState, m_meshObj,
			   m_burgers2d_icPulse,
			   m_burgers2d_icSpread,
			   m_burgers2d_x0,
			   m_burgers2d_y0);
      }

    return initialState;
  }

protected:
  void initializeJacobian(jacobian_type & J)
  {
    J.resize(m_numDofSampleMesh, m_numDofStencilMesh);

    using Tr = Eigen::Triplet<scalar_type>;
    std::vector<Tr> trList;

    initializeJacobianForInnerCells(trList);
    if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersDirichlet){
      initializeJacobianForNearBoundaryCells(trList);
    }

    J.setFromTriplets(trList.begin(), trList.end());

    // compress to make it Csr
    if (!J.isCompressed()){
      J.makeCompressed();
    }
  }

  // note that we MUST template U_t, V_t because
  // when doing bindings, these are deduced to be a Eigen Ref
  template<class U_t, class V_t>
  void velocityAndOptionalJacobian(const U_t & U,
				   const scalar_type currentTime,
				   V_t & V,
				   jacobian_type * J) const
  {

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp parallel
    {
#endif

      // for omp, these are private variables for each thread
      // edge reconstructions
      edge_rec_type uMinusHalfNeg(m_numDofPerCell);
      edge_rec_type uMinusHalfPos(m_numDofPerCell);
      edge_rec_type uPlusHalfNeg(m_numDofPerCell);
      edge_rec_type uPlusHalfPos(m_numDofPerCell);
      // fluxes
      flux_type fluxL(m_numDofPerCell);
      flux_type fluxF(m_numDofPerCell);;
      flux_type fluxR(m_numDofPerCell);
      flux_type fluxB(m_numDofPerCell);

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
      ::pressiodemoapps::set_zero_omp(V);
#else
      ::pressiodemoapps::set_zero(V);
#endif

      if (J){
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
	::pressiodemoapps::set_zero_omp(*J);
#else
	::pressiodemoapps::set_zero(*J);
#endif
      }

      fillGhosts(U);

      if (J){
	if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersDirichlet){
	  velocityAndJacobianDirichletImpl(U, currentTime, V, *J,
					   fluxL, fluxF, fluxR, fluxB,
					   uMinusHalfNeg, uMinusHalfPos,
					   uPlusHalfNeg,  uPlusHalfPos);
	}
	else{
	  velocityAndJacobianPeriodicImpl(U, currentTime, V, *J,
					   fluxL, fluxF, fluxR, fluxB,
					   uMinusHalfNeg, uMinusHalfPos,
					   uPlusHalfNeg,  uPlusHalfPos);
	}
      }
      else{

	if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersDirichlet){
	  velocityOnlyNearBdCellsImpl(U, currentTime, V,
				      fluxL, fluxF, fluxR, fluxB,
				      uMinusHalfNeg, uMinusHalfPos,
				      uPlusHalfNeg,  uPlusHalfPos);
	}

	velocityOnlyInnerCellsImpl(U, currentTime, V,
				   fluxL, fluxF, fluxR, fluxB,
				   uMinusHalfNeg, uMinusHalfPos,
				   uPlusHalfNeg,  uPlusHalfPos);
      }

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
}//end omp parallel
#endif

  }

private:
  template<class Tr>
  void initializeJacobianForInnerCells(std::vector<Tr> & trList)
  {
    // for inner cells, the Jacobian is exact and depends on
    // the scheme wanted by the user, no special treatment needed

    const scalar_type zero = 0;
    const auto & graph = m_meshObj.graph();
    // only grab the graph rows for INNER cells (i.e. AWAY from boundaries)
    const auto & targetGraphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (std::size_t it=0; it<targetGraphRows.size(); ++it)
      {
	const auto smPt = targetGraphRows[it];

	// find out which row in the jacobian we are dealing with
	const auto jacRowOfCurrCellFirstDof = smPt*m_numDofPerCell;

	// initialize jacobian block entries wrt current cell's dofs
	const auto jacColOfCurrCellRho = graph(smPt, 0)*m_numDofPerCell;
	for (int k=0; k<m_numDofPerCell; ++k){
	  for (int j=0; j<m_numDofPerCell; ++j){
	    trList.push_back( Tr(jacRowOfCurrCellFirstDof+k,
				 jacColOfCurrCellRho+j,
				 zero) );
	  }
	}

	// now fill wrt neighbors: this depends on the scheme
	// so I need to use the largest between the stencil needed
	// for inviscid and viscous reconstructions.
	// Most likely the inviscid always wins so use that for now.
	// if assert breaks, this will need to be fixed.
	assert(reconstructionTypeToStencilSize(m_inviscidFluxRecEn)
	       >= reconstructionTypeToStencilSize(m_viscousFluxRecEn));

	// find out how many neighboring cells: recall this is 2d
	const int numNeighbors =
	  (m_inviscidFluxRecEn == InviscidFluxReconstruction::FirstOrder) ? 4
	    : (m_inviscidFluxRecEn == InviscidFluxReconstruction::Weno3) ? 8
	      : (m_inviscidFluxRecEn == InviscidFluxReconstruction::Weno5) ? 12 : -1;
	assert(numNeighbors > 0);

	for (int i=1; i<=numNeighbors; ++i){
	  const auto colInd = graph(smPt, i)*m_numDofPerCell;
	  for (int k=0; k<m_numDofPerCell; ++k){
	    for (int j=0; j<m_numDofPerCell; ++j){
	      trList.push_back( Tr(jacRowOfCurrCellFirstDof+k, colInd+j, zero) );
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
	const auto jacRowOfCurrCellFirstDof = smPt*m_numDofPerCell;
	const auto jacColOfCurrCellFirstDof = graph(smPt, 0)*m_numDofPerCell;

	// wrt current cell's dofs
	for (int k=0; k<m_numDofPerCell; ++k){
	  for (int j=0; j<m_numDofPerCell; ++j){
	    trList.push_back( Tr(jacRowOfCurrCellFirstDof+k, jacColOfCurrCellFirstDof+j, zero) );
	  }
	}

	// wrt neighbors
	// for near-bd, we only do first-order Jacobian for now
	for (int i=1; i<=4; ++i){
	  const auto nID = graph(smPt, i);
	  if (nID != -1){
	    const auto ci = nID*m_numDofPerCell;
	    for (int k=0; k<m_numDofPerCell; ++k){
	      for (int j=0; j<m_numDofPerCell; ++j){
		trList.push_back( Tr(jacRowOfCurrCellFirstDof+k, ci+j, zero) );
	      }
	    }
	  }
	}
      }
  }

  template<class U_t>
  void fillGhosts(const U_t & U) const
  {
    const auto stencilSize = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersDirichlet)
    {
      using ghost_filler_t  = ::pressiodemoapps::impladvdiff2d::BurgersDiriFiller<
	U_t, MeshType, ghost_container_type>;
      ghost_filler_t ghF(stencilSize, U, m_meshObj,
			 m_ghostLeft, m_ghostFront,
			 m_ghostRight, m_ghostBack);

      const auto & rowsBd = m_meshObj.graphRowsOfCellsNearBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
      for (int it=0; it<rowsBd.size(); ++it){
	ghF(rowsBd[it], it);
      }
    }

    else{
      // no op
    }
  }

  template<class U_t, class V_t>
  void velocityAndJacobianDirichletImpl(const U_t & U,
					const scalar_type currentTime,
					V_t & V,
					jacobian_type & J,
					flux_type & fluxL,
					flux_type & fluxF,
					flux_type & fluxR,
					flux_type & fluxB,
					edge_rec_type & uMinusHalfNeg,
					edge_rec_type & uMinusHalfPos,
					edge_rec_type & uPlusHalfNeg,
					edge_rec_type & uPlusHalfPos) const
  {

    // flux jacobians
    flux_jac_type fluxJacLNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacLPos(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacFNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacFPos(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacRNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacRPos(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacBNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacBPos(m_numDofPerCell, m_numDofPerCell);

    int nonZerosCountBeforeComputing = 0;
    nonZerosCountBeforeComputing = J.nonZeros();

    // near boundary I have be careful because
    // the jacobian can only be first order for now
    if (m_inviscidFluxRecEn == InviscidFluxReconstruction::FirstOrder){
      velocityAndJacNearBDCellsImplFirstOrder(U, currentTime, V, J,
					      fluxL, fluxF, fluxR, fluxB,
					      fluxJacLNeg, fluxJacLPos,
					      fluxJacFNeg, fluxJacFPos,
					      fluxJacRNeg, fluxJacRPos,
					      fluxJacBNeg, fluxJacBPos,
					      uMinusHalfNeg, uMinusHalfPos,
					      uPlusHalfNeg,  uPlusHalfPos);
    }
    else{
      velocityAndJacNearBDCellsImplDifferentScheme(U, currentTime, V, J,
						   fluxL, fluxF, fluxR, fluxB,
						   fluxJacLNeg, fluxJacLPos,
						   fluxJacFNeg, fluxJacFPos,
						   fluxJacRNeg, fluxJacRPos,
						   fluxJacBNeg, fluxJacBPos,
						   uMinusHalfNeg, uMinusHalfPos,
						   uPlusHalfNeg,  uPlusHalfPos);
    }

    velocityAndJacInnerCellsImpl(U, currentTime, V, J,
				 fluxL, fluxF, fluxR, fluxB,
				 fluxJacLNeg, fluxJacLPos,
				 fluxJacFNeg, fluxJacFPos,
				 fluxJacRNeg, fluxJacRPos,
				 fluxJacBNeg, fluxJacBPos,
				 uMinusHalfNeg, uMinusHalfPos,
				 uPlusHalfNeg,  uPlusHalfPos);

    assert(J.nonZeros() == nonZerosCountBeforeComputing);
  }

  template<class U_t, class V_t>
  void velocityAndJacobianPeriodicImpl(const U_t & U,
				       const scalar_type currentTime,
				       V_t & V,
				       jacobian_type & J,
				       flux_type & fluxL,
				       flux_type & fluxF,
				       flux_type & fluxR,
				       flux_type & fluxB,
				       edge_rec_type & uMinusHalfNeg,
				       edge_rec_type & uMinusHalfPos,
				       edge_rec_type & uPlusHalfNeg,
				       edge_rec_type & uPlusHalfPos) const
  {

    // flux jacobians
    flux_jac_type fluxJacLNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacLPos(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacFNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacFPos(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacRNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacRPos(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacBNeg(m_numDofPerCell, m_numDofPerCell);
    flux_jac_type fluxJacBPos(m_numDofPerCell, m_numDofPerCell);

    int nonZerosCountBeforeComputing = 0;
    nonZerosCountBeforeComputing = J.nonZeros();

    // this is for periodic so we don't have boundary cells,
    // all cells count as "inner cells"
    velocityAndJacInnerCellsImpl(U, currentTime, V, J,
				 fluxL, fluxF, fluxR, fluxB,
				 fluxJacLNeg, fluxJacLPos,
				 fluxJacFNeg, fluxJacFPos,
				 fluxJacRNeg, fluxJacRPos,
				 fluxJacBNeg, fluxJacBPos,
				 uMinusHalfNeg, uMinusHalfPos,
				 uPlusHalfNeg,  uPlusHalfPos);

    assert(J.nonZeros() == nonZerosCountBeforeComputing);
  }

  template<class U_t, class V_t>
  void velocityAndJacNearBDCellsImplFirstOrder(const U_t & U,
					       const scalar_type currentTime,
					       V_t & V,
					       jacobian_type & J,
					       flux_type & fluxL,
					       flux_type & fluxF,
					       flux_type & fluxR,
					       flux_type & fluxB,
					       flux_jac_type & fluxJacLNeg,
					       flux_jac_type & fluxJacLPos,
					       flux_jac_type & fluxJacFNeg,
					       flux_jac_type & fluxJacFPos,
					       flux_jac_type & fluxJacRNeg,
					       flux_jac_type & fluxJacRPos,
					       flux_jac_type & fluxJacBNeg,
					       flux_jac_type & fluxJacBPos,
					       edge_rec_type & uMinusHalfNeg,
					       edge_rec_type & uMinusHalfPos,
					       edge_rec_type & uPlusHalfNeg,
					       edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    assert(m_inviscidFluxRecEn == InviscidFluxReconstruction::FirstOrder);

    // if here, then the scheme for velocity matches
    // the one for Jacobian so we can use same functors
    const auto stencilSize = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    stencil_container_type stencilVals(m_numDofPerCell*stencilSize);

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, stencil_container_type, U_t, MeshType, ghost_container_type>;
    stencil_filler_t FillStencilX(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				   U, m_meshObj, m_ghostLeft, m_ghostRight,
				   stencilVals, xAxis);

    stencil_filler_t FillStencilY(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				   U, m_meshObj, m_ghostBack, m_ghostFront,
				   stencilVals, yAxis);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impl::ComputeDirectionalFluxBalanceFirstOrderJacobianOnBoundaryCell<
	  pda::impladvdiff2d::ComputeDirectionalFluxValuesAndJacobians<
	    pda::impl::ReconstructorFromStencil<
	      edge_rec_type, stencil_container_type>,
	    scalar_type, flux_type, flux_jac_type>,
	  dimensionality, MeshType, jacobian_type>,
      V_t, scalar_type
      >;

    functor_type funcx(V, m_meshObj.dxInv(),
		       /* end args for velo */
		       J, xAxis, m_meshObj,
		       /* end args for jac */
		       m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR,
		       fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
		       /* end args for flux */
		       toReconstructionScheme(m_inviscidFluxRecEn), stencilVals,
		       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		       /* end args for reconstructor */
		       );

    functor_type funcy(V, m_meshObj.dyInv(),
		       /* end args for velo */
		       J, yAxis, m_meshObj,
		       /* end args for jac */
		       m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF,
		       fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
		       /* end args for flux */
		       toReconstructionScheme(m_inviscidFluxRecEn), stencilVals,
		       uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		       /* end args for reconstructor */
		       );

    std::array<scalar_type, 2> bdCellJacFactorsX;
    std::array<scalar_type, 2> bdCellJacFactorsY;
    bdCellJacFactorsX.fill(static_cast<scalar_type>(-1));
    bdCellJacFactorsY.fill(static_cast<scalar_type>(-1));

    const auto & graph      = m_meshObj.graph();
    const auto & graphRows = m_meshObj.graphRowsOfCellsNearBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];

      FillStencilX(smPt, it, m_numDofPerCell);
      funcx(smPt, m_numDofPerCell, bdCellJacFactorsX, 2);
      FillStencilY(smPt, it, m_numDofPerCell);
      funcy(smPt, m_numDofPerCell, bdCellJacFactorsY, 2);

      if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersDirichlet)
	{
	  constexpr auto two  = static_cast<scalar_type>(2);
	  const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
	  const auto dyInvSq  = m_meshObj.dyInv()*m_meshObj.dyInv();
	  const auto diffDxInvSq  = m_burgers2d_diffusion*dxInvSq;
	  const auto diffDyInvSq  = m_burgers2d_diffusion*dyInvSq;

	  const auto cellIndexLeft  = graph(smPt, 1);
	  const auto cellIndexFront = graph(smPt, 2);
	  const auto cellIndexRight = graph(smPt, 3);
	  const auto cellIndexBack  = graph(smPt, 4);

	  const auto vIndex = smPt*m_numDofPerCell;
	  const auto uIndex = graph(smPt, 0)*m_numDofPerCell;

	  auto selfValue1 = -two*diffDxInvSq -two*diffDyInvSq;
	  auto selfValue2 = -two*diffDxInvSq -two*diffDyInvSq;

	  if (cellIndexLeft != -1){
	    const auto uIndexLeft  = cellIndexLeft*m_numDofPerCell;
	    J.coeffRef(vIndex, uIndexLeft) += diffDxInvSq;
	    J.coeffRef(vIndex+1, uIndexLeft+1) += diffDxInvSq;
	  }else{
	    selfValue1 += -diffDxInvSq;
	    selfValue2 += -diffDxInvSq;
	  }

	  if (cellIndexFront != -1){
	  const auto uIndexFront = cellIndexFront*m_numDofPerCell;
	    J.coeffRef(vIndex, uIndexFront) += diffDyInvSq;
	    J.coeffRef(vIndex+1, uIndexFront+1) += diffDyInvSq;
	  }else{
	    selfValue1 += -diffDyInvSq;
	    selfValue2 += -diffDyInvSq;
	  }

	  if (cellIndexRight != -1){
	    const auto uIndexRight = cellIndexRight*m_numDofPerCell;
	    J.coeffRef(vIndex, uIndexRight) += diffDxInvSq;
	    J.coeffRef(vIndex+1, uIndexRight+1) += diffDxInvSq;
	  }
	  else{
	    selfValue1 += -diffDxInvSq;
	    selfValue2 += -diffDxInvSq;
	  }

	  if (cellIndexBack != -1){
	  const auto uIndexBack  = cellIndexBack*m_numDofPerCell;
	    J.coeffRef(vIndex, uIndexBack) += diffDyInvSq;
	    J.coeffRef(vIndex+1, uIndexBack+1) += diffDyInvSq;
	  }
	  else{
	    selfValue1 += -diffDyInvSq;
	    selfValue2 += -diffDyInvSq;
	  }

	  J.coeffRef(vIndex, uIndex) += selfValue1;
	  J.coeffRef(vIndex+1, uIndex+1) += selfValue2;
      }
    }
  }

  template<class U_t, class V_t>
  void velocityAndJacNearBDCellsImplDifferentScheme(const U_t & U,
						    const scalar_type currentTime,
						    V_t & V,
						    jacobian_type & J,
						    flux_type & fluxL,
						    flux_type & fluxF,
						    flux_type & fluxR,
						    flux_type & fluxB,
						    flux_jac_type & fluxJacLNeg,
						    flux_jac_type & fluxJacLPos,
						    flux_jac_type & fluxJacFNeg,
						    flux_jac_type & fluxJacFPos,
						    flux_jac_type & fluxJacRNeg,
						    flux_jac_type & fluxJacRPos,
						    flux_jac_type & fluxJacBNeg,
						    flux_jac_type & fluxJacBPos,
						    edge_rec_type & uMinusHalfNeg,
						    edge_rec_type & uMinusHalfPos,
						    edge_rec_type & uPlusHalfNeg,
						    edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;

    // if here, then the velocity must be computed with Weno,
    /// while the jacobian must be computed with first order

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, stencil_container_type,
      U_t, MeshType, ghost_container_type>;

    // *****************************
    // *** functors for velocity ***
    // *****************************
    const auto stencilSizeForV = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    stencil_container_type stencilValsForV(m_numDofPerCell*stencilSizeForV);

    stencil_filler_t FillStencilVeloX(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				      U, m_meshObj, m_ghostLeft, m_ghostRight,
				      stencilValsForV, xAxis);
    stencil_filler_t FillStencilVeloY(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				      U, m_meshObj, m_ghostBack, m_ghostFront,
				      stencilValsForV, yAxis);

    using velo_functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impladvdiff2d::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  scalar_type, flux_type>,
      V_t, scalar_type
      >;

    velo_functor_type funcVeloX(V, m_meshObj.dxInv(),
				/* end args for velo */
				m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR,
				/* end args for flux */
				toReconstructionScheme(m_inviscidFluxRecEn), stencilValsForV,
				uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
				/* end args for reconstructor */
				);

    velo_functor_type funcVeloY(V, m_meshObj.dyInv(),
				/* end args for velo */
				m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF,
				/* end args for flux */
				toReconstructionScheme(m_inviscidFluxRecEn), stencilValsForV,
				uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
				/* end args for reconstructor */
				);

    // *****************************
    // *** functors for jacobian ***
    // *****************************
    const auto firstOrderRec = pda::InviscidFluxReconstruction::FirstOrder;
    const auto stencilSizeForJ = reconstructionTypeToStencilSize(firstOrderRec);
    stencil_container_type stencilValsForJ(m_numDofPerCell*stencilSizeForJ);
    stencil_filler_t FillStencilJacX(stencilSizeForJ,
				     U, m_meshObj, m_ghostLeft, m_ghostRight,
				     stencilValsForJ, xAxis);
    stencil_filler_t FillStencilJacY(stencilSizeForJ,
				     U, m_meshObj, m_ghostBack, m_ghostFront,
				     stencilValsForJ, yAxis);

    using jac_functor_type =
      pda::impl::ComputeDirectionalFluxBalanceFirstOrderJacobianOnBoundaryCell<
	pda::impladvdiff2d::ComputeDirectionalFluxJacobians<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  scalar_type, flux_jac_type>,
      dimensionality, MeshType, jacobian_type
      >;

    jac_functor_type funcJacX(J, xAxis, m_meshObj,
			      /* end args for jac */
			      m_probEn, m_inviscidFluxSchemeEn, normalX_,
			      fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
			      /* end args for flux */
			      toReconstructionScheme(firstOrderRec), stencilValsForJ,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			      /* end args for reconstructor */
			      );

    jac_functor_type funcJacY(J, yAxis, m_meshObj,
			      /* end args for jac */
			      m_probEn, m_inviscidFluxSchemeEn, normalY_,
			      fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
			      /* end args for flux */
			      toReconstructionScheme(firstOrderRec), stencilValsForJ,
			      uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
			      /* end args for reconstructor */
			      );

    std::array<scalar_type, 2> bdCellJacFactorsX;
    std::array<scalar_type, 2> bdCellJacFactorsY;
    bdCellJacFactorsX.fill(static_cast<scalar_type>(-1));
    bdCellJacFactorsY.fill(static_cast<scalar_type>(-1));

    // ************
    // loop
    // ************
    const auto & graph      = m_meshObj.graph();
    const auto & graphRows = m_meshObj.graphRowsOfCellsNearBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {

      const auto smPt = graphRows[it];
      const auto vIndex = smPt*m_numDofPerCell;

      FillStencilVeloX(smPt, it, m_numDofPerCell);
      funcVeloX(smPt, m_numDofPerCell);
      FillStencilJacX(smPt, it, m_numDofPerCell);
      funcJacX(smPt, m_numDofPerCell, bdCellJacFactorsX, 1);
      if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersDirichlet){
	constexpr auto two  = static_cast<scalar_type>(2);
	const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
	const auto diffDxInvSq  = m_burgers2d_diffusion*dxInvSq;

	if ( stencilSizeForV == 3){
	  V(vIndex)   += diffDxInvSq*( stencilValsForV(4)
				       - two*stencilValsForV(2) + stencilValsForV(0) );
	  V(vIndex+1) += diffDxInvSq*( stencilValsForV(5)
				       - two*stencilValsForV(3) + stencilValsForV(1) );
	}
	else if ( stencilSizeForV == 5){
	  V(vIndex)   += diffDxInvSq*( stencilValsForV(6)
				       - two*stencilValsForV(4) + stencilValsForV(2) );
	  V(vIndex+1) += diffDxInvSq*( stencilValsForV(7)
				       - two*stencilValsForV(5) + stencilValsForV(3) );
	}
	else if ( stencilSizeForV == 7){
	  V(vIndex)   += diffDxInvSq*( stencilValsForV(8)
				       - two*stencilValsForV(6) + stencilValsForV(4) );
	  V(vIndex+1) += diffDxInvSq*( stencilValsForV(9)
				       - two*stencilValsForV(7) + stencilValsForV(5) );
	}
      }

      FillStencilVeloY(smPt, it, m_numDofPerCell);
      funcVeloY(smPt, m_numDofPerCell);
      FillStencilJacY(smPt, it, m_numDofPerCell);
      funcJacY(smPt, m_numDofPerCell, bdCellJacFactorsY, 1);
      if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersDirichlet){
	constexpr auto two  = static_cast<scalar_type>(2);
	const auto dyInvSq  = m_meshObj.dyInv()*m_meshObj.dyInv();
	const auto diffDyInvSq  = m_burgers2d_diffusion*dyInvSq;

	if ( stencilSizeForV == 3){
	  V(vIndex)   += diffDyInvSq*( stencilValsForV(4)
				       - two*stencilValsForV(2) + stencilValsForV(0) );
	  V(vIndex+1) += diffDyInvSq*( stencilValsForV(5)
				       - two*stencilValsForV(3) + stencilValsForV(1) );
	}
	else if ( stencilSizeForV == 5){
	  V(vIndex)   += diffDyInvSq*( stencilValsForV(6)
				       - two*stencilValsForV(4) + stencilValsForV(2) );
	  V(vIndex+1) += diffDyInvSq*( stencilValsForV(7)
				       - two*stencilValsForV(5) + stencilValsForV(3) );
	}
	else if ( stencilSizeForV == 7){
	  V(vIndex)   += diffDyInvSq*( stencilValsForV(8)
				       - two*stencilValsForV(6) + stencilValsForV(4) );
	  V(vIndex+1) += diffDyInvSq*( stencilValsForV(9)
				       - two*stencilValsForV(7) + stencilValsForV(5) );
	}
      }

      if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersDirichlet)
	{
	  constexpr auto two  = static_cast<scalar_type>(2);
	  const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
	  const auto dyInvSq  = m_meshObj.dyInv()*m_meshObj.dyInv();
	  const auto diffDxInvSq  = m_burgers2d_diffusion*dxInvSq;
	  const auto diffDyInvSq  = m_burgers2d_diffusion*dyInvSq;

	  const auto cellIndexLeft  = graph(smPt, 1);
	  const auto cellIndexFront = graph(smPt, 2);
	  const auto cellIndexRight = graph(smPt, 3);
	  const auto cellIndexBack  = graph(smPt, 4);

	  const auto vIndex = smPt*m_numDofPerCell;
	  const auto uIndex = graph(smPt, 0)*m_numDofPerCell;

	  auto selfValue1 = -two*diffDxInvSq -two*diffDyInvSq;
	  auto selfValue2 = -two*diffDxInvSq -two*diffDyInvSq;

	  if (cellIndexLeft != -1){
	    const auto uIndexLeft  = graph(smPt, 1)*m_numDofPerCell;
	    J.coeffRef(vIndex, uIndexLeft) += diffDxInvSq;
	    J.coeffRef(vIndex+1, uIndexLeft+1) += diffDxInvSq;
	  }else{
	    selfValue1 += -diffDxInvSq;
	    selfValue2 += -diffDxInvSq;
	  }

	  if (cellIndexFront != -1){
	    const auto uIndexFront = graph(smPt, 2)*m_numDofPerCell;
	    J.coeffRef(vIndex, uIndexFront) += diffDyInvSq;
	    J.coeffRef(vIndex+1, uIndexFront+1) += diffDyInvSq;
	  }else{
	    selfValue1 += -diffDyInvSq;
	    selfValue2 += -diffDyInvSq;
	  }

	  if (cellIndexRight != -1){
	    const auto uIndexRight = graph(smPt, 3)*m_numDofPerCell;
	    J.coeffRef(vIndex, uIndexRight) += diffDxInvSq;
	    J.coeffRef(vIndex+1, uIndexRight+1) += diffDxInvSq;
	  }
	  else{
	    selfValue1 += -diffDxInvSq;
	    selfValue2 += -diffDxInvSq;
	  }

	  if (cellIndexBack != -1){
	    const auto uIndexBack  = graph(smPt, 4)*m_numDofPerCell;
	    J.coeffRef(vIndex, uIndexBack) += diffDyInvSq;
	    J.coeffRef(vIndex+1, uIndexBack+1) += diffDyInvSq;
	  }
	  else{
	    selfValue1 += -diffDyInvSq;
	    selfValue2 += -diffDyInvSq;
	  }

	  J.coeffRef(vIndex, uIndex) += selfValue1;
	  J.coeffRef(vIndex+1, uIndex+1) += selfValue2;
	}
    }
  }


  template<class U_t, class V_t>
  void velocityAndJacInnerCellsImpl(const U_t & U,
				    const scalar_type currentTime,
				    V_t & V,
				    jacobian_type & J,
				    flux_type & fluxL,
				    flux_type & fluxF,
				    flux_type & fluxR,
				    flux_type & fluxB,
				    flux_jac_type & fluxJacLNeg,
				    flux_jac_type & fluxJacLPos,
				    flux_jac_type & fluxJacFNeg,
				    flux_jac_type & fluxJacFPos,
				    flux_jac_type & fluxJacRNeg,
				    flux_jac_type & fluxJacRPos,
				    flux_jac_type & fluxJacBNeg,
				    flux_jac_type & fluxJacBPos,
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

    // allocate gradients of reconstructed states
    // the size depends on the scheme selected
    const auto stencilSize = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    assert(stencilSize >= reconstructionTypeToStencilSize(m_viscousFluxRecEn));

    reconstruction_gradient_t gradLNeg(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradLPos(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradFNeg(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradFPos(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRNeg(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradRPos(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradBNeg(m_numDofPerCell, stencilSize-1);
    reconstruction_gradient_t gradBPos(m_numDofPerCell, stencilSize-1);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impl::ComputeDirectionalFluxBalanceJacobianOnInteriorCell<
	  pda::impladvdiff2d::ComputeDirectionalFluxValuesAndJacobians<
	    pda::impl::ReconstructorForDiscreteFunction<
	      dimensionality, MeshType, U_t, edge_rec_type, reconstruction_gradient_t>,
	    scalar_type, flux_type, flux_jac_type>,
	  dimensionality, MeshType, jacobian_type>,
      V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.dxInv(),
		    /* end args for velo */
		    J, xAxis, m_meshObj,
		    /* end args for jac */
		    m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR,
		    fluxJacLNeg, fluxJacLPos, fluxJacRNeg, fluxJacRPos,
		    /* end args for flux */
		    xAxis, toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		    gradLNeg, gradLPos, gradRNeg, gradRPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.dyInv(),
		    /* end args for velo */
		    J, yAxis, m_meshObj,
		    /* end args for jac */
		    m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF,
		    fluxJacBNeg, fluxJacBPos, fluxJacFNeg, fluxJacFPos,
		    /* end args for flux */
		    yAxis, toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos,
		    gradBNeg, gradBPos, gradFNeg, gradFPos
		    /* end args for reconstructor */
		    );

    const auto & graphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];
      Fx(smPt, m_numDofPerCell);
      Fy(smPt, m_numDofPerCell);

      if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersPeriodic ||
	  m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersDirichlet){
	addBurgersDiffusionAndSourceToVelocityAndJacobianInnerCells(U, V, J, smPt);
      }
    }
  }

  template<class U_t, class V_t>
  void velocityOnlyNearBdCellsImpl(const U_t & U,
				   const scalar_type currentTime,
				   V_t & V,
				   flux_type & fluxL,
				   flux_type & fluxF,
				   flux_type & fluxR,
				   flux_type & fluxB,
				   edge_rec_type & uMinusHalfNeg,
				   edge_rec_type & uMinusHalfPos,
				   edge_rec_type & uPlusHalfNeg,
				   edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;

    const auto stencilSize = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    assert(stencilSize >= reconstructionTypeToStencilSize(m_viscousFluxRecEn));
    stencil_container_type stencilVals(m_numDofPerCell*stencilSize);

    using stencil_filler_t  = pda::impl::StencilFiller<
      dimensionality, stencil_container_type, U_t, MeshType, ghost_container_type>;
    stencil_filler_t FillStencilX(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				  U, m_meshObj, m_ghostLeft, m_ghostRight,
				  stencilVals, xAxis);

    stencil_filler_t FillStencilY(reconstructionTypeToStencilSize(m_inviscidFluxRecEn),
				  U, m_meshObj, m_ghostBack, m_ghostFront,
				  stencilVals, yAxis);

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impladvdiff2d::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorFromStencil<
	    edge_rec_type, stencil_container_type>,
	  scalar_type, flux_type>,
      V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.dxInv(),
		    /* end args for velo */
		    m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR,
		    /* end args for flux */
		    toReconstructionScheme(m_inviscidFluxRecEn), stencilVals,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.dyInv(),
		    /* end args for velo */
		    m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF,
		    /* end args for flux */
		    toReconstructionScheme(m_inviscidFluxRecEn), stencilVals,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg, uPlusHalfPos
		    /* end args for reconstructor */
		    );

    const auto & graphRows = m_meshObj.graphRowsOfCellsNearBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it)
    {
      const auto smPt = graphRows[it];
      const auto vIndex = smPt*m_numDofPerCell;

      FillStencilX(smPt, it, m_numDofPerCell);
      Fx(smPt, m_numDofPerCell);
      if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersDirichlet){
	constexpr auto two  = static_cast<scalar_type>(2);
	const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
	const auto diffDxInvSq  = m_burgers2d_diffusion*dxInvSq;
	V(vIndex)   += diffDxInvSq*( stencilVals(4) - two*stencilVals(2) + stencilVals(0) );
	V(vIndex+1) += diffDxInvSq*( stencilVals(5) - two*stencilVals(3) + stencilVals(1) );
      }

      FillStencilY(smPt, it, m_numDofPerCell);
      Fy(smPt, m_numDofPerCell);
      if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersDirichlet){
	constexpr auto two  = static_cast<scalar_type>(2);
	const auto dyInvSq  = m_meshObj.dyInv()*m_meshObj.dyInv();
	const auto diffDyInvSq  = m_burgers2d_diffusion*dyInvSq;
	V(vIndex)   += diffDyInvSq*( stencilVals(4) - two*stencilVals(2) + stencilVals(0) );
	V(vIndex+1) += diffDyInvSq*( stencilVals(5) - two*stencilVals(3) + stencilVals(1) );
      }

    }
  }

  template<class U_t, class V_t>
  void velocityOnlyInnerCellsImpl(const U_t & U,
				  const scalar_type currentTime,
				  V_t & V,
				  flux_type & fluxL,
				  flux_type & fluxF,
				  flux_type & fluxR,
				  flux_type & fluxB,
				  edge_rec_type & uMinusHalfNeg,
				  edge_rec_type & uMinusHalfPos,
				  edge_rec_type & uPlusHalfNeg,
				  edge_rec_type & uPlusHalfPos) const
  {
    namespace pda = ::pressiodemoapps;
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;

    using functor_type =
      pda::impl::ComputeDirectionalFluxBalance<
	pda::impladvdiff2d::ComputeDirectionalFluxValues<
	  pda::impl::ReconstructorForDiscreteFunction<
	    dimensionality, MeshType, U_t, edge_rec_type>,
	  scalar_type, flux_type>,
      V_t, scalar_type
      >;

    functor_type Fx(V, m_meshObj.dxInv(),
		    /* end args for velo */
		    m_probEn, m_inviscidFluxSchemeEn, normalX_, fluxL, fluxR,
		    /* end args for flux */
		    xAxis, toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		    /* end args for reconstructor */
		    );

    functor_type Fy(V, m_meshObj.dyInv(),
		    /* end args for velo */
		    m_probEn, m_inviscidFluxSchemeEn, normalY_, fluxB, fluxF,
		    /* end args for flux */
		    yAxis, toReconstructionScheme(m_inviscidFluxRecEn), U, m_meshObj,
		    uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
		    /* end args for reconstructor */
		    );

    const auto & graphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
    for (decltype(graphRows.size()) it=0; it<graphRows.size(); ++it){
      const auto smPt = graphRows[it];
      Fx(smPt, m_numDofPerCell);
      Fy(smPt, m_numDofPerCell);

      if (m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersPeriodic ||
	  m_probEn == ::pressiodemoapps::AdvectionDiffusion2d::BurgersDirichlet){
	addBurgersDiffusionToVelocityInnerCells(U, V, smPt);
      }
    }
  }

  template<class U_t, class V_t>
  void addBurgersDiffusionToVelocityInnerCells(const U_t & U,
					       V_t & V,
					       index_t smPt) const
  {
    constexpr auto two  = static_cast<scalar_type>(2);
    const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto dyInvSq  = m_meshObj.dyInv()*m_meshObj.dyInv();
    const auto diffDxInvSq  = m_burgers2d_diffusion*dxInvSq;
    const auto diffDyInvSq  = m_burgers2d_diffusion*dyInvSq;

    const auto & graph     = m_meshObj.graph();
    const auto vIndex      = smPt*m_numDofPerCell;
    const auto uIndex      = graph(smPt, 0)*m_numDofPerCell;
    const auto uIndexLeft  = graph(smPt, 1)*m_numDofPerCell;
    const auto uIndexFront = graph(smPt, 2)*m_numDofPerCell;
    const auto uIndexRight = graph(smPt, 3)*m_numDofPerCell;
    const auto uIndexBack  = graph(smPt, 4)*m_numDofPerCell;
    V(vIndex) += diffDxInvSq*( U(uIndexRight) - two*U(uIndex) + U(uIndexLeft) );
    V(vIndex) += diffDyInvSq*( U(uIndexFront) - two*U(uIndex) + U(uIndexBack) );

    V(vIndex+1) += diffDxInvSq*( U(uIndexRight+1) - two*U(uIndex+1) + U(uIndexLeft+1) );
    V(vIndex+1) += diffDyInvSq*( U(uIndexFront+1) - two*U(uIndex+1) + U(uIndexBack+1) );
  }

  template<class U_t, class V_t>
  void addBurgersDiffusionAndSourceToVelocityAndJacobianInnerCells(const U_t & U,
								   V_t & V,
								   jacobian_type & J,
								   index_t smPt) const
  {
    constexpr auto two  = static_cast<scalar_type>(2);
    const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto dyInvSq  = m_meshObj.dyInv()*m_meshObj.dyInv();
    const auto diffDxInvSq  = m_burgers2d_diffusion*dxInvSq;
    const auto diffDyInvSq  = m_burgers2d_diffusion*dyInvSq;

    const auto & graph     = m_meshObj.graph();
    const auto vIndex      = smPt*m_numDofPerCell;
    const auto uIndex      = graph(smPt, 0)*m_numDofPerCell;
    const auto uIndexLeft  = graph(smPt, 1)*m_numDofPerCell;
    const auto uIndexFront = graph(smPt, 2)*m_numDofPerCell;
    const auto uIndexRight = graph(smPt, 3)*m_numDofPerCell;
    const auto uIndexBack  = graph(smPt, 4)*m_numDofPerCell;
    V(vIndex) += diffDxInvSq*( U(uIndexRight) - two*U(uIndex) + U(uIndexLeft) );
    V(vIndex) += diffDyInvSq*( U(uIndexFront) - two*U(uIndex) + U(uIndexBack) );

    V(vIndex+1) += diffDxInvSq*( U(uIndexRight+1) - two*U(uIndex+1) + U(uIndexLeft+1) );
    V(vIndex+1) += diffDyInvSq*( U(uIndexFront+1) - two*U(uIndex+1) + U(uIndexBack+1) );

    J.coeffRef(vIndex, uIndex) += -two*diffDxInvSq - two*diffDyInvSq;
    J.coeffRef(vIndex, uIndexLeft)  += diffDxInvSq;
    J.coeffRef(vIndex, uIndexFront) += diffDyInvSq;
    J.coeffRef(vIndex, uIndexRight) += diffDxInvSq;
    J.coeffRef(vIndex, uIndexBack)  += diffDyInvSq;
    J.coeffRef(vIndex+1, uIndex+1) += -two*diffDxInvSq - two*diffDyInvSq;
    J.coeffRef(vIndex+1, uIndexLeft+1)  += diffDxInvSq;
    J.coeffRef(vIndex+1, uIndexFront+1) += diffDyInvSq;
    J.coeffRef(vIndex+1, uIndexRight+1) += diffDxInvSq;
    J.coeffRef(vIndex+1, uIndexBack+1)  += diffDyInvSq;
  }

private:
  void allocateGhosts()
  {
    const auto stencilSize    = reconstructionTypeToStencilSize(m_inviscidFluxRecEn);
    const auto numGhostValues = m_numDofPerCell*((stencilSize-1)/2);

    const index_t s1 = m_meshObj.numCellsBd();
    ::pressiodemoapps::resize(m_ghostLeft, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostFront,s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight,s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostBack, s1, numGhostValues);
  }

protected:
  // common to all problems
  int m_numDofPerCell = {};
  ::pressiodemoapps::AdvectionDiffusion2d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_inviscidFluxRecEn;
  ::pressiodemoapps::InviscidFluxScheme m_inviscidFluxSchemeEn;
  ::pressiodemoapps::ViscousFluxReconstruction m_viscousFluxRecEn;

  const MeshType & m_meshObj;
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  mutable ghost_container_type m_ghostLeft;
  mutable ghost_container_type m_ghostFront;
  mutable ghost_container_type m_ghostRight;
  mutable ghost_container_type m_ghostBack;

  const std::array<scalar_type, 2> normalX_{1, 0};
  const std::array<scalar_type, 2> normalY_{0, 1};

  // parameters specific to problems
  // will need to handle this better later
  scalar_type m_burgers2d_icPulse = {};
  scalar_type m_burgers2d_icSpread = {};
  scalar_type m_burgers2d_diffusion = {};
  scalar_type m_burgers2d_x0 = {};
  scalar_type m_burgers2d_y0 = {};
};

template<class MeshType> constexpr int EigenApp<MeshType>::dimensionality;

}}//end namespace
#endif
