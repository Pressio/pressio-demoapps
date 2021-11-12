
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

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
#include "Eigen/Sparse"
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#include <omp.h>
#endif

namespace pressiodemoapps{ namespace ee{ namespace impl{

template<
  class ScalarType,
  class MeshType,
  class StateType,
  class VelocityType,
  class GhostContainerType
  >
class Euler2dAppRhsOnly
{

public:
  using index_t	       = typename MeshType::index_t;
  using scalar_type    = ScalarType;
  using state_type     = StateType;
  using velocity_type  = VelocityType;

  using ghost_vals_container_type   = GhostContainerType;
  using stencil_vals_container_type = state_type;
  using flux_t		 = state_type;
  using edge_rec_t	 = state_type;

  static constexpr int     dimensionality{2};
  static constexpr index_t numDofPerCell{4};

public:
  Euler2dAppRhsOnly(const MeshType & meshObj,
		    ::pressiodemoapps::Euler2d probEn,
		    ::pressiodemoapps::InviscidFluxReconstruction recEn,
		    ::pressiodemoapps::InviscidFluxScheme fluxEnum,
		    int icIdentifier)
    : m_meshObj(meshObj)
    ,m_probEn(probEn)
    ,m_recEn(recEn)
    ,m_fluxEn(fluxEnum)
    ,m_icIdentifier(icIdentifier)
#if defined PRESSIODEMOAPPS_ENABLE_BINDINGS
    ,m_stencilVals(1)
    // note that when doing bindings, I need to first construct
    // ghost with {1,1} just so that the numpy array picks up they
    // are 2dim array otherwise it thinks they are 1d array.
    // The right allocation for these is then done inside allocateGhosts.
    ,m_ghostLeft({1,1})
    ,m_ghostFront({1,1})
    ,m_ghostRight({1,1})
    ,m_ghostBack({1,1})
#endif
  {
    computeDofs();
    allocateStencilValuesContainer();
    allocateGhosts();
  }

  state_type initialCondition() const{
    return initialConditionImpl();
  }

  scalar_type gamma()           const{ return m_gamma; }
  index_t totalDofSampleMesh()  const{ return m_numDofSampleMesh; }
  index_t totalDofStencilMesh() const{ return m_numDofStencilMesh; }

  velocity_type createVelocity() const{
    velocity_type V(m_numDofSampleMesh);
    return V;
  }

  void velocity(const state_type & state,
		const scalar_type currentTime,
		velocity_type & V) const
  {
    flux_t FLeft(numDofPerCell);
    flux_t FRight(numDofPerCell);
    flux_t FFront(numDofPerCell);
    flux_t FBack(numDofPerCell);
    edge_rec_t uMinusHalfNeg(numDofPerCell);
    edge_rec_t uMinusHalfPos(numDofPerCell);
    edge_rec_t uPlusHalfNeg (numDofPerCell);
    edge_rec_t uPlusHalfPos (numDofPerCell);

    fillGhosts(state, currentTime);
    velocityCellsNearBdImpl(state, currentTime, V,
			    FLeft, FRight, FBack, FFront,
			    uMinusHalfNeg, uMinusHalfPos,
			    uPlusHalfNeg,  uPlusHalfPos);
    velocityInnerCellsImpl(state, currentTime, V,
			   FLeft, FRight, FBack, FFront,
			   uMinusHalfNeg, uMinusHalfPos,
			   uPlusHalfNeg,  uPlusHalfPos);
  }

private:
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

protected:
  void velocityInnerCellsImpl(const state_type & U,
			      const scalar_type currentTime,
			      velocity_type & V,
			      flux_t & FLeft,
			      flux_t & FRight,
			      flux_t & FBack,
			      flux_t & FFront,
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
      dimensionality, edge_rec_t, state_type, MeshType>;

    rec_fnct_t ReconstructorX(xAxis, m_recEn, U, m_meshObj,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);

    rec_fnct_t ReconstructorY(yAxis, m_recEn, U, m_meshObj,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);

    const auto & rowsIn = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (int it=0; it<rowsIn.size(); ++it)
    {
      const auto smPt = rowsIn[it];

      // X
      ReconstructorX.template operator()<numDofPerCell>(smPt);
      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFourDof(FLeft, uMinusHalfNeg, uMinusHalfPos, normalX_, m_gamma);
	  eeRusanovFluxFourDof(FRight, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gamma);
	  break;
	}

      // Y
      ReconstructorY.template operator()<numDofPerCell>(smPt);
      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFourDof(FBack, uMinusHalfNeg, uMinusHalfPos, normalY_, m_gamma);
	  eeRusanovFluxFourDof(FFront, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gamma);
	  break;
	}

      const auto vIndex = smPt*numDofPerCell;
      V(vIndex)   = dxInv*(FLeft(0) - FRight(0)) + dyInv*(FBack(0) - FFront(0));
      V(vIndex+1) = dxInv*(FLeft(1) - FRight(1)) + dyInv*(FBack(1) - FFront(1));
      V(vIndex+2) = dxInv*(FLeft(2) - FRight(2)) + dyInv*(FBack(2) - FFront(2));
      V(vIndex+3) = dxInv*(FLeft(3) - FRight(3)) + dyInv*(FBack(3) - FFront(3));
    }
  }

  void velocityCellsNearBdImpl(const state_type & U,
			       const scalar_type currentTime,
			       velocity_type & V,
			       flux_t & FLeft,
			       flux_t & FRight,
			       flux_t & FBack,
			       flux_t & FFront,
			       edge_rec_t & uMinusHalfNeg,
			       edge_rec_t & uMinusHalfPos,
			       edge_rec_t & uPlusHalfNeg,
			       edge_rec_t & uPlusHalfPos) const
  {

    constexpr int xAxis = 1;
    constexpr int yAxis = 2;

    using sfiller_t = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_vals_container_type, state_type, MeshType, ghost_vals_container_type>;
    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromStencilFourDofPerCell
      <edge_rec_t, stencil_vals_container_type>;

    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    sfiller_t StencilFillerX(stencilSize, U, m_meshObj,
			     m_ghostLeft, m_ghostRight,
			     m_stencilVals, xAxis);

    sfiller_t StencilFillerY(stencilSize, U, m_meshObj,
			     m_ghostBack, m_ghostFront,
			     m_stencilVals, yAxis);

    rec_fnct_t Reconstructor(m_recEn, m_stencilVals,
			     uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg,  uPlusHalfPos);

    const auto dxInv   = m_meshObj.dxInv();
    const auto dyInv   = m_meshObj.dyInv();
    const auto & specialRows = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<specialRows.size(); ++it)
    {
      const auto smPt    = specialRows[it];

      // X
      StencilFillerX(smPt, it);
      Reconstructor();
      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFourDof(FLeft, uMinusHalfNeg, uMinusHalfPos, normalX_, m_gamma);
	  eeRusanovFluxFourDof(FRight, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gamma);
	  break;
	}

      // Y
      StencilFillerY(smPt, it);
      Reconstructor();
      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFourDof(FBack, uMinusHalfNeg, uMinusHalfPos, normalY_, m_gamma);
	  eeRusanovFluxFourDof(FFront, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gamma);
	  break;
	}

      const auto vIndex = smPt*numDofPerCell;
      V(vIndex)   = dxInv*(FLeft(0) - FRight(0)) + dyInv*(FBack(0) - FFront(0));
      V(vIndex+1) = dxInv*(FLeft(1) - FRight(1)) + dyInv*(FBack(1) - FFront(1));
      V(vIndex+2) = dxInv*(FLeft(2) - FRight(2)) + dyInv*(FBack(2) - FFront(2));
      V(vIndex+3) = dxInv*(FLeft(3) - FRight(3)) + dyInv*(FBack(3) - FFront(3));
    }
  }

  void fillGhosts(const state_type & U, const scalar_type currentTime) const
  {
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    if (m_probEn == ::pressiodemoapps::Euler2d::SedovFull or
	m_probEn == ::pressiodemoapps::Euler2d::Riemann or
	m_probEn == ::pressiodemoapps::Euler2d::testingonlyneumann)
    {
      using ghost_filler_t  = ::pressiodemoapps::impl::Ghost2dNeumannFiller<
	numDofPerCell, state_type, MeshType, ghost_vals_container_type>;
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
	state_type, MeshType, ghost_vals_container_type>;
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
	state_type, MeshType, ghost_vals_container_type>;
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
	state_type, MeshType, ghost_vals_container_type>;
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

  mutable stencil_vals_container_type m_stencilVals;
  mutable ghost_vals_container_type m_ghostLeft;
  mutable ghost_vals_container_type m_ghostFront;
  mutable ghost_vals_container_type m_ghostRight;
  mutable ghost_vals_container_type m_ghostBack;

  const std::array<scalar_type, 2> normalX_{1, 0};
  const std::array<scalar_type, 2> normalY_{0, 1};
};



#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class ScalarType, class MeshType>
using EigenEuler2dAppRhsOnly = Euler2dAppRhsOnly<
  ScalarType, MeshType,
  Eigen::Matrix<ScalarType,-1,1>,
  Eigen::Matrix<ScalarType,-1,1>,
  Eigen::Matrix<ScalarType,-1,-1, Eigen::RowMajor>
  >;

template<class ScalarType, class MeshType>
class EigenEuler2dAppWithJacobian
  : public EigenEuler2dAppRhsOnly<ScalarType, MeshType>
{
  using base_t = EigenEuler2dAppRhsOnly<ScalarType, MeshType>;
  using base_t::numDofPerCell;

public:
  using typename base_t::index_t;
  using typename base_t::scalar_type;
  using typename base_t::state_type;
  using typename base_t::velocity_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, int>;

  using typename base_t::ghost_vals_container_type;
  using typename base_t::stencil_vals_container_type;
  using typename base_t::flux_t;
  using typename base_t::edge_rec_t;
  using flux_jac_type = Eigen::Matrix<ScalarType, numDofPerCell, numDofPerCell>;

public:
  template<class ...Args>
  EigenEuler2dAppWithJacobian(Args && ... args)
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

  void velocity(const state_type & state,
		const scalar_type currentTime,
		velocity_type & V) const
  {

    flux_jac_type JLneg, JLpos;
    flux_jac_type JRneg, JRpos;
    flux_jac_type JFneg, JFpos;
    flux_jac_type JBneg, JBpos;

    // flux at faces
    flux_t FLeft(numDofPerCell);
    flux_t FRight(numDofPerCell);
    flux_t FFront(numDofPerCell);
    flux_t FBack(numDofPerCell);

    // reconstructed values at faces
    edge_rec_t uMinusHalfNeg(numDofPerCell);
    edge_rec_t uMinusHalfPos(numDofPerCell);
    edge_rec_t uPlusHalfNeg(numDofPerCell);
    edge_rec_t uPlusHalfPos(numDofPerCell);

    base_t::fillGhosts(state, currentTime);

    if (m_onlyComputeVelocity){
      base_t::velocityCellsNearBdImpl(state, currentTime, V,
				      FLeft, FRight, FBack, FFront,
				      uMinusHalfNeg, uMinusHalfPos,
				      uPlusHalfNeg,  uPlusHalfPos);
      base_t::velocityInnerCellsImpl(state, currentTime, V,
				     FLeft, FRight, FBack, FFront,
				     uMinusHalfNeg, uMinusHalfPos,
				     uPlusHalfNeg,  uPlusHalfPos);
    }

    else{

      auto values = m_jacobian.valuePtr();
      for (int i=0; i<m_jacobian.nonZeros(); ++i){
	values[i] = 0;
      }

      velocityAndJacInnerCellsImpl(state, currentTime, V,
				   FLeft, FRight, FBack, FFront,
				   JLneg, JLpos,
				   JRneg, JRpos,
				   JBneg, JBpos,
				   JFneg, JFpos,
				   uMinusHalfNeg, uMinusHalfPos,
				   uPlusHalfNeg,  uPlusHalfPos);

      velocityAndJacNearBDCellsImpl(state, currentTime, V,
				   FLeft, FRight, FBack, FFront,
				   JLneg, JLpos,
				   JRneg, JRpos,
				   JBneg, JBpos,
				   JFneg, JFpos,
				   uMinusHalfNeg, uMinusHalfPos,
				   uPlusHalfNeg,  uPlusHalfPos);


#ifdef NDEBUG
      // ensure that nonzeros count does not change
      // this can happen for instance if during the Jacobian
      // evlauation, new non-zero elements get inserted
      assert(m_jacNonZerosCount == m_jacobian.nonZeros());
#endif
    }
  }

  jacobian_type createJacobian() const {
    return m_jacobian;
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
    using Tr = Eigen::Triplet<ScalarType>;
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

    m_jacobian.setFromTriplets(trList.begin(), trList.end());
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

  void velocityAndJacInnerCellsImpl(const state_type & U,
				    const scalar_type currentTime,
				    velocity_type & V,
				    flux_t & FLeft,
				    flux_t & FRight,
				    flux_t & FBack,
				    flux_t & FFront,
				    flux_jac_type & JLneg,
				    flux_jac_type & JLpos,
				    flux_jac_type & JRneg,
				    flux_jac_type & JRpos,
				    flux_jac_type & JBneg,
				    flux_jac_type & JBpos,
				    flux_jac_type & JFneg,
				    flux_jac_type & JFpos,
				    edge_rec_t & uMinusHalfNeg,
				    edge_rec_t & uMinusHalfPos,
				    edge_rec_t & uPlusHalfNeg,
				    edge_rec_t & uPlusHalfPos) const
  {
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    const auto dxInv   = m_meshObj.dxInv();
    const auto dyInv   = m_meshObj.dyInv();

    using reconstructor_functor_t = ::pressiodemoapps::impl::ReconstructorFromState<
      base_t::dimensionality, edge_rec_t, state_type, MeshType>;

    using jac_fnct_t = ::pressiodemoapps::impl::FirstOrderInnerCellJacobianFunctor<
      base_t::dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, MeshType>;

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
    const auto firstOrder = ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder;
    edge_rec_t uMinusHalfNegForJ(numDofPerCell);
    edge_rec_t uMinusHalfPosForJ(numDofPerCell);
    edge_rec_t uPlusHalfNegForJ(numDofPerCell);
    edge_rec_t uPlusHalfPosForJ(numDofPerCell);
    reconstructor_functor_t ReconstructorXForJ(xAxis, firstOrder, U, m_meshObj,
					       uMinusHalfNegForJ, uMinusHalfPosForJ,
					       uPlusHalfNegForJ,  uPlusHalfPosForJ);

    reconstructor_functor_t ReconstructorYForJ(yAxis, firstOrder, U, m_meshObj,
					       uMinusHalfNegForJ, uMinusHalfPosForJ,
					       uPlusHalfNegForJ,  uPlusHalfPosForJ);

    jac_fnct_t CellJacobianFunctorX(m_jacobian, m_meshObj, JLneg, JLpos, JRneg, JRpos, xAxis);
    jac_fnct_t CellJacobianFunctorY(m_jacobian, m_meshObj, JBneg, JBpos, JFneg, JFpos, yAxis);

    const auto & targetGraphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (int it=0; it<targetGraphRows.size(); ++it)
    {
      const auto smPt = targetGraphRows[it];
      const auto vIndex = smPt*numDofPerCell;

      // *** X ***
      ReconstructorX.template operator()<numDofPerCell>(smPt);

      // note that REGARDLESS of the reconstruction scheme,
      // we currently only have only first-order Jacobian so we need
      // to run the reconstructor for the Jacobian
      // which will ensure that uMinusNegForJ, etc have the right values
      ReconstructorXForJ.template operator()<numDofPerCell>(smPt);

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFourDof(FLeft,  uMinusHalfNeg, uMinusHalfPos, normalX_, m_gamma);
	  eeRusanovFluxFourDof(FRight, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gamma);

	  eeRusanovFluxJacobianFourDof(JLneg, JLpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
				       normalX_, m_gamma);
	  eeRusanovFluxJacobianFourDof(JRneg, JRpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
				       normalX_, m_gamma);
	  break;
	}

      CellJacobianFunctorX(smPt);

      // *** Y ***
      ReconstructorY.template operator()<numDofPerCell>(smPt);

      // note that REGARDLESS of the reconstruction scheme,
      // we currently only have only first-order Jacobian so we need
      // to run the reconstructor for the Jacobian
      // which will ensure that uMinusNegForJ, etc have the right values
      ReconstructorYForJ.template operator()<numDofPerCell>(smPt);

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFourDof(FBack,  uMinusHalfNeg, uMinusHalfPos, normalY_, m_gamma);
	  eeRusanovFluxFourDof(FFront, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gamma);

	  eeRusanovFluxJacobianFourDof(JBneg, JBpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
				       normalY_, m_gamma);
	  eeRusanovFluxJacobianFourDof(JFneg, JFpos, uPlusHalfNegForJ,  uPlusHalfPosForJ,
				       normalY_, m_gamma);
	  break;
	}

      CellJacobianFunctorY(smPt);

      V(vIndex)   = dxInv*(FLeft(0) - FRight(0)) + dyInv*(FBack(0) - FFront(0));
      V(vIndex+1) = dxInv*(FLeft(1) - FRight(1)) + dyInv*(FBack(1) - FFront(1));
      V(vIndex+2) = dxInv*(FLeft(2) - FRight(2)) + dyInv*(FBack(2) - FFront(2));
      V(vIndex+3) = dxInv*(FLeft(3) - FRight(3)) + dyInv*(FBack(3) - FFront(3));
    }
  }

  void velocityAndJacNearBDCellsImpl(const state_type & U,
				     const scalar_type currentTime,
				     velocity_type & V,
				     flux_t & FLeft,
				     flux_t & FRight,
				     flux_t & FBack,
				     flux_t & FFront,
				     flux_jac_type & JLneg,
				     flux_jac_type & JLpos,
				     flux_jac_type & JRneg,
				     flux_jac_type & JRpos,
				     flux_jac_type & JBneg,
				     flux_jac_type & JBpos,
				     flux_jac_type & JFneg,
				     flux_jac_type & JFpos,
				     edge_rec_t & uMinusHalfNeg,
				     edge_rec_t & uMinusHalfPos,
				     edge_rec_t & uPlusHalfNeg,
				     edge_rec_t & uPlusHalfPos) const
  {

    constexpr int xAxis = 1;
    constexpr int yAxis = 2;
    const auto dxInv   = m_meshObj.dxInv();
    const auto dyInv   = m_meshObj.dyInv();

    using stencil_filler_t = ::pressiodemoapps::impl::StencilFiller<
      base_t::dimensionality, numDofPerCell, stencil_vals_container_type,
      state_type, MeshType, ghost_vals_container_type>;
    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromStencilFourDofPerCell
      <edge_rec_t, stencil_vals_container_type>;

    using jac_fnct_t = ::pressiodemoapps::impl::FirstOrderBdCellJacobianFunctor<
      base_t::dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, MeshType>;

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
    const auto cellJacOrdEn = ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder;

    stencil_vals_container_type stencilValsForJ = {};
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

    std::array<scalar_type, numDofPerCell> bcCellJacFactorsDefault;
    std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveX;
    std::array<scalar_type, numDofPerCell> bcCellJacFactorsReflectiveY;
    bcCellJacFactorsDefault.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsReflectiveX.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsReflectiveY.fill(static_cast<scalar_type>(1));
    bcCellJacFactorsReflectiveX[1] = static_cast<ScalarType>(-1);
    bcCellJacFactorsReflectiveY[2] = static_cast<ScalarType>(-1);

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

      FillStencilValuesXFunctorForJ(smPt, it);
      FaceValuesReconstructFunctorForJ();

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFourDof(FLeft, uMinusHalfNeg, uMinusHalfPos, normalX_, m_gamma);
	  eeRusanovFluxFourDof(FRight, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gamma);

	  eeRusanovFluxJacobianFourDof(JLneg, JLpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
				       normalX_, m_gamma);
	  eeRusanovFluxJacobianFourDof(JRneg, JRpos, uPlusHalfNegForJ, uPlusHalfPosForJ,
				       normalX_, m_gamma);
	  break;
	}

      auto bcType = findCellBdType(smPt, xAxis);
      const auto & factorsX = (bcType == 1) ? bcCellJacFactorsReflectiveX : bcCellJacFactorsDefault;
      CellJacobianFunctorX(smPt, factorsX, bcType);

      // *** Y ***
      StencilFillerY(smPt, it);
      Reconstructor();

      FillStencilValuesYFunctorForJ(smPt, it);
      FaceValuesReconstructFunctorForJ();

      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFourDof(FBack, uMinusHalfNeg, uMinusHalfPos, normalY_, m_gamma);
	  eeRusanovFluxFourDof(FFront, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gamma);

	  eeRusanovFluxJacobianFourDof(JBneg, JBpos, uMinusHalfNegForJ, uMinusHalfPosForJ,
				       normalY_, m_gamma);
	  eeRusanovFluxJacobianFourDof(JFneg, JFpos, uPlusHalfNegForJ,  uPlusHalfPosForJ,
				       normalY_, m_gamma);
	  break;
	}

      bcType = findCellBdType(smPt, yAxis);
      const auto & factorsY = (bcType == 1) ? bcCellJacFactorsReflectiveY : bcCellJacFactorsDefault;
      CellJacobianFunctorY(smPt, factorsY, bcType);

      V(vIndex)   = dxInv*(FLeft(0) - FRight(0)) + dyInv*(FBack(0) - FFront(0));
      V(vIndex+1) = dxInv*(FLeft(1) - FRight(1)) + dyInv*(FBack(1) - FFront(1));
      V(vIndex+2) = dxInv*(FLeft(2) - FRight(2)) + dyInv*(FBack(2) - FFront(2));
      V(vIndex+3) = dxInv*(FLeft(3) - FRight(3)) + dyInv*(FBack(3) - FFront(3));
    }
  }


  using base_t::m_meshObj;
  using base_t::m_ghostLeft;
  using base_t::m_ghostFront;
  using base_t::m_ghostRight;
  using base_t::m_ghostBack;
  using base_t::m_stencilVals;
  using base_t::normalX_;
  using base_t::normalY_;
  using base_t::m_gamma;
  using base_t::m_numDofStencilMesh;
  using base_t::m_numDofSampleMesh;
  using base_t::m_recEn;
  using base_t::m_probEn;
  using base_t::m_fluxEn;
  mutable jacobian_type m_jacobian = {};
  std::size_t m_jacNonZerosCount = {};
  bool m_onlyComputeVelocity = false;
};
#endif

}}}//end namespace
#endif
