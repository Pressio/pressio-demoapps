
#ifndef PRESSIODEMOAPPS_EULER1D_APP_HPP_
#define PRESSIODEMOAPPS_EULER1D_APP_HPP_

// note that the code below is pretty ugly for now, but it works.
// this is inside impl namespace for a reason, and will need to be
// improved later on but we have a starting point.

#include "../eulerCommon/fluxes.hpp"
#include "../eulerCommon/jacobians.hpp"
#include "../eulerCommon/rankine_hugoniot.hpp"
#include "./initial_condition.hpp"
#include "./ghost_filler.hpp"
#include "../stencil_filler.hpp"
#include "../reconstructor_from_stencil.hpp"
#include "../reconstructor_from_state.hpp"

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
#include "Eigen/Sparse"
#endif

namespace pressiodemoapps{ namespace impl{

template<
  class ScalarType,
  class MeshType,
  class StateType,
  class VelocityType,
  class GhostContainerType
  >
class Euler1dAppRhsOnly
{

public:
  using index_t	       = typename MeshType::index_t;
  using scalar_type    = ScalarType;
  using state_type     = StateType;
  using velocity_type  = VelocityType;

  using ghost_vals_container_type   = GhostContainerType;
  using stencil_vals_container_type = state_type;
  using flux_t			= state_type;
  using edge_rec_t		= state_type;

  static constexpr int dimensionality{1};
  static constexpr index_t numDofPerCell{3};

  // // stencil filler functor type
  // using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
  //   dimensionality, numDofPerCell, stencil_vals_container_type,
  //   state_type, MeshType, ghost_vals_container_type>;

  // // reconstructor functor type
  // using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromStencilThreeDofPerCell<
  //   edge_rec_t, stencil_vals_container_type>;

public:
  Euler1dAppRhsOnly(const MeshType & meshObj,
	      ::pressiodemoapps::Euler1d probEnum,
	      ::pressiodemoapps::InviscidFluxReconstruction recEnum,
	      ::pressiodemoapps::InviscidFluxScheme fluxEnum)
    : m_meshObj(meshObj), m_probEn(probEnum), m_recEn(recEnum), m_fluxEn(fluxEnum)
  {

    // calculate total num of dofs on sample and stencil mesh
    m_numDofStencilMesh = m_meshObj.stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * numDofPerCell;

    allocateStencilValuesContainer();
    allocateGhostValuesContainer();
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

  scalar_type gamma()		const{ return m_gamma; }
  index_t totalDofSampleMesh()  const{ return m_numDofSampleMesh; }
  index_t totalDofStencilMesh() const{ return m_numDofStencilMesh; }

  velocity_type createVelocity() const {
    velocity_type V(m_numDofSampleMesh);
    return V;
  }

  void velocity(const state_type & state,
		const scalar_type currentTime,
		velocity_type & V) const
  {
    flux_t FL(numDofPerCell);
    flux_t FR(numDofPerCell);
    edge_rec_t uMinusHalfNeg(numDofPerCell);
    edge_rec_t uMinusHalfPos(numDofPerCell);
    edge_rec_t uPlusHalfNeg (numDofPerCell);
    edge_rec_t uPlusHalfPos (numDofPerCell);

    fillGhostsIfNeeded(state, currentTime);
    velocityCellsNearBdImpl(state, currentTime,
			    V, FL, FR,
			    uMinusHalfNeg, uMinusHalfPos,
			    uPlusHalfNeg,  uPlusHalfPos);

    velocityInnerCellsImpl(state, currentTime,
			   V, FL, FR,
			   uMinusHalfNeg, uMinusHalfPos,
			   uPlusHalfNeg,  uPlusHalfPos);
  }

private:
  void fillGhostsIfNeeded(const state_type & U,
			  scalar_type /*currTime*/) const
  {
    const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_recEn);
    // only need ghosts for specific problems
    if (m_probEn == ::pressiodemoapps::Euler1d::Sod or
	m_probEn == ::pressiodemoapps::Euler1d::Lax)
    {
      using ghost_filler_t  = ::pressiodemoapps::impl::Ghost1dNeumannFiller<
	numDofPerCell, state_type, MeshType, ghost_vals_container_type>;
      ghost_filler_t ghF(stencilSizeNeeded, U, m_meshObj, m_ghostLeft, m_ghostRight);
      ghF();
    }
  }

  void velocityInnerCellsImpl(const state_type & U,
			      const scalar_type currentTime,
			      velocity_type & V,
			      flux_t & FL,
			      flux_t & FR,
			      edge_rec_t & uMinusHalfNeg,
			      edge_rec_t & uMinusHalfPos,
			      edge_rec_t & uPlusHalfNeg,
			      edge_rec_t & uPlusHalfPos) const
  {
    /*
      for inner cells, we do not need to worry aboit boundaries,
      so we can index the state directly to to reconstruction
      and therefore we do not need to fill a temporary stencil array values
    */

    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromState<
      dimensionality, edge_rec_t, state_type, MeshType>;

    rec_fnct_t Reconstructor(m_recEn, U, m_meshObj,
			     uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg,  uPlusHalfPos);

    // deal with cells away from boundaries
    const auto dxInv = m_meshObj.dxInv();
    const auto & rowsIn = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (int it=0; it<rowsIn.size(); ++it)
    {
      const auto smPt = rowsIn[it];

      Reconstructor.template operator()<numDofPerCell>(smPt);
      switch(m_fluxEn)
      {
      case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  ee::impl::eeRusanovFluxThreeDof(FL, uMinusHalfNeg, uMinusHalfPos, m_gamma);
	  ee::impl::eeRusanovFluxThreeDof(FR, uPlusHalfNeg,  uPlusHalfPos,  m_gamma);
	  break;
      }

      const auto vIndex = smPt*numDofPerCell;
      V(vIndex)   = dxInv*(FL(0) - FR(0));
      V(vIndex+1) = dxInv*(FL(1) - FR(1));
      V(vIndex+2) = dxInv*(FL(2) - FR(2));
    }
  }

  void velocityCellsNearBdImpl(const state_type & U,
			       const scalar_type currentTime,
			       velocity_type & V,
			       flux_t & FL,
			       flux_t & FR,
			       edge_rec_t & uMinusHalfNeg,
			       edge_rec_t & uMinusHalfPos,
			       edge_rec_t & uPlusHalfNeg,
			       edge_rec_t & uPlusHalfPos) const
  {
    /*
      for cells that are near boundaries, we need use the stencil
      values to handle properly the ghost cells.
      Once we fill the stencil values needed, we can do
      reconstruction and fluxes.
     */

    using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_vals_container_type,
      state_type, MeshType, ghost_vals_container_type>;
    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromStencilThreeDofPerCell<
      edge_rec_t, stencil_vals_container_type>;

    const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_recEn);
    sfiller_t StencilFiller(stencilSizeNeeded, U, m_meshObj,
			    m_ghostLeft, m_ghostRight, m_stencilVals);
    rec_fnct_t Reconstructor(m_recEn, m_stencilVals,
			     uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg,  uPlusHalfPos);

    const auto dxInv = m_meshObj.dxInv();
    const auto & specialRows = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<specialRows.size(); ++it)
    {
      const auto smPt = specialRows[it];

      StencilFiller(smPt);
      Reconstructor();
      switch(m_fluxEn)
      {
      case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  ee::impl::eeRusanovFluxThreeDof(FL, uMinusHalfNeg, uMinusHalfPos, m_gamma);
	  ee::impl::eeRusanovFluxThreeDof(FR, uPlusHalfNeg,  uPlusHalfPos,  m_gamma);
	  break;
      }

      const auto vIndex = smPt*numDofPerCell;
      V(vIndex)   = dxInv*(FL(0) - FR(0));
      V(vIndex+1) = dxInv*(FL(1) - FR(1));
      V(vIndex+2) = dxInv*(FL(2) - FR(2));
    }
  }

private:
  void allocateStencilValuesContainer()
  {
    // the stencil size needed is determined by the desired reconstruction
    // kind NOT from the mesh. THis is important because for example
    // the mesh can have a wider connectivity that what is needed.
    const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_recEn);
    ::pressiodemoapps::resize(m_stencilVals, numDofPerCell*stencilSizeNeeded);
  }

  void allocateGhostValuesContainer()
  {
    /*
      In 1d, the ghost values at left and right are stored in arrays.
      The size of this array depends on the stencil.

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

    const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_recEn);
    const auto ghostStorageSize = numDofPerCell*((stencilSizeNeeded-1)/2);
    ::pressiodemoapps::resize(m_ghostLeft,  ghostStorageSize);
    ::pressiodemoapps::resize(m_ghostRight, ghostStorageSize);
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

  mutable stencil_vals_container_type m_stencilVals = {};
  mutable ghost_vals_container_type   m_ghostLeft = {};
  mutable ghost_vals_container_type   m_ghostRight = {};
};


#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
template<class ScalarType, class MeshType>
using EigenEuler1dAppRhsOnly = Euler1dAppRhsOnly<ScalarType, MeshType,
						 Eigen::Matrix<ScalarType,-1,1>,
						 Eigen::Matrix<ScalarType,-1,1>,
						 Eigen::Matrix<ScalarType,-1,1>>;

template<class ScalarType, class MeshType>
class EigenEuler1dAppTWithJacobian
  : public Euler1dAppRhsOnly<ScalarType, MeshType,
			      Eigen::Matrix<ScalarType,-1,1>,
			      Eigen::Matrix<ScalarType,-1,1>,
			      Eigen::Matrix<ScalarType,-1,1>>
{
  using base_t = Euler1dAppRhsOnly<ScalarType, MeshType,
				    Eigen::Matrix<ScalarType,-1,1>,
				    Eigen::Matrix<ScalarType,-1,1>,
				    Eigen::Matrix<ScalarType,-1,1>>;

  using base_t::numDofPerCell;

public:
  using typename base_t::index_t;
  using typename base_t::scalar_type;
  using typename base_t::state_type;
  using typename base_t::velocity_type;
  using jacobian_type  = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, int>;

private:
  using typename base_t::ghost_vals_container_type;
  using typename base_t::stencil_vals_container_type;
  using typename base_t::flux_t;
  using typename base_t::edge_rec_t;
  using flux_jac_type = Eigen::Matrix<ScalarType,3,3>;
  using Tr = Eigen::Triplet<ScalarType>;

public:
  template<class ...Args>
  EigenEuler1dAppTWithJacobian(Args && ... args)
    : base_t(std::forward<Args>(args)...)
  {}

  void fuseVelocityAndJacobian() {
    m_fuseVelocityAndJacobian = true;
  }

  void velocity(const state_type & state,
		const scalar_type currentTime,
		velocity_type & V) const
  {
    base_t::fillGhostsIfNeeded(state, currentTime);

    // jac at left face wrt to left,right recontructes states
    flux_jac_type JLneg, JLpos;
    // jac at right face wrt to left,right recontructes states
    flux_jac_type JRneg, JRpos;

    flux_t FL(numDofPerCell);
    flux_t FR(numDofPerCell);
    edge_rec_t uMinusHalfNeg(numDofPerCell);
    edge_rec_t uMinusHalfPos(numDofPerCell);
    edge_rec_t uPlusHalfNeg (numDofPerCell);
    edge_rec_t uPlusHalfPos (numDofPerCell);

    // handle all inner cells
    velocityAndJacInnerCellsImpl(state, currentTime,
				 V, FL, FR, JLneg, JLpos, JRneg, JRpos,
				 uMinusHalfNeg, uMinusHalfPos,
				 uPlusHalfNeg,  uPlusHalfPos);

    // handle all near-BD cells
    velocityAndJacCellsNearBdImpl(state, currentTime,
				  V, FL, FR, JLneg, JLpos, JRneg, JRpos,
				  uMinusHalfNeg, uMinusHalfPos,
				  uPlusHalfNeg,  uPlusHalfPos);

    // handle all BD cells
    // functor to compute the cell Jacobian for 1st,weno3 (this is what
    // we have done with Patrick on phone)
  }

  void jacobian(const state_type & state,
		const scalar_type timeValue,
		jacobian_type & J) const
  {
    // assuming jacobian has been computed above
    J = m_jacobian;
  }

private:

  void velocityAndJacInnerCellsImpl(const state_type & U,
				    const scalar_type currentTime,
				    velocity_type & V,
				    flux_t & FL,
				    flux_t & FR,
				    flux_jac_type & JLneg,
				    flux_jac_type & JLpos,
				    flux_jac_type & JRneg,
				    flux_jac_type & JRpos,
				    edge_rec_t & uMinusHalfNeg,
				    edge_rec_t & uMinusHalfPos,
				    edge_rec_t & uPlusHalfNeg,
				    edge_rec_t & uPlusHalfPos) const
  {

  }

  void velocityAndJacNearBDCellsImpl(const state_type & U,
				     const scalar_type currentTime,
				     velocity_type & V,
				     flux_t & FL,
				     flux_t & FR,
				     flux_jac_type & JLneg,
				     flux_jac_type & JLpos,
				     flux_jac_type & JRneg,
				     flux_jac_type & JRpos,
				     edge_rec_t & uMinusHalfNeg,
				     edge_rec_t & uMinusHalfPos,
				     edge_rec_t & uPlusHalfNeg,
				     edge_rec_t & uPlusHalfPos) const
  {

    using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_vals_container_type,
      state_type, MeshType, ghost_vals_container_type>;
    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromStencilThreeDofPerCell<
      edge_rec_t, stencil_vals_container_type>;

    const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_recEn);
    sfiller_t FillStencilValuesFunctor(stencilSizeNeeded, U, m_meshObj,
				       m_ghostLeft, m_ghostRight,
				       m_stencilVals);
    rec_fnct_t FaceValuesReconstructFunctor(base_t::m_recEn, m_stencilVals,
					    uMinusHalfNeg, uMinusHalfPos,
					    uPlusHalfNeg,  uPlusHalfPos);

    using jac_fnct_t = FirstOrderNearBDCellJacobianThreeDofFunctor<
      dimensionality, scalar_type, jacobian_type, MeshType>;

    jac_fnct_t CellJacobianFunctor(m_meshObj, jacobian, JLneg, JLpos, LRneg, JRpos);

    const auto dxInv = m_meshObj.dxInv();
    const auto & specialRows = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<specialRows.size(); ++it)
    {
      const auto smPt = specialRows[it];

      StencilFiller(smPt);
      Reconstructor();
      switch(m_fluxEn)
      {
      case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	ee::impl::eeRusanovFluxThreeDof(FL, uMinusHalfNeg, uMinusHalfPos, m_gamma);
	ee::impl::eeRusanovFluxThreeDof(FR, uPlusHalfNeg,  uPlusHalfPos,  m_gamma);

	ee::impl::eeRusanovFluxJacobianThreeDof(JLneg, JLpos, uMinusHalfNeg, uMinusHalfPos, m_gamma);
	ee::impl::eeRusanovFluxJacobianThreeDof(JRneg, JRpos, uPlusHalfNeg,  uPlusHalfPos, m_gamma);
	break;
      }

      CellJacobianFunctor(smPt);

      const auto vIndex = smPt*numDofPerCell;
      V(vIndex)   = dxInv*(FL(0) - FR(0));
      V(vIndex+1) = dxInv*(FL(1) - FR(1));
      V(vIndex+2) = dxInv*(FL(2) - FR(2));
    }
  }

  void velocityAndJacobianFirstOrderImplInnnerCells(const state_type & U,
						    const scalar_type timeValue,
						    velocity_type & V) const
  {
    // // ----------
    // // 1. prepare
    // // ----------
    // const auto stencilSizeNeeded = reconstructionTypeToStencilSize(base_t::m_recEn);

    // // jac at left face wrt to left,right recontructes states
    // flux_jac_type JLneg, JLpos;
    // // jac at right face wrt to left,right recontructes states
    // flux_jac_type JRneg, JRpos;

    // flux_t FL(numDofPerCell);
    // flux_t FR(numDofPerCell);
    // edge_rec_t uMinusHalfNeg(numDofPerCell);
    // edge_rec_t uMinusHalfPos(numDofPerCell);
    // edge_rec_t uPlusHalfNeg (numDofPerCell);
    // edge_rec_t uPlusHalfPos (numDofPerCell);

    // // -------------------
    // // 4. create functors
    // // -------------------
    // typename base_t::sfiller_t FillStencilValuesFunctor(stencilSizeNeeded, U, m_meshObj,
    // 							m_ghostLeft, m_ghostRight,
    // 							m_stencilVals);
    // typename base_t::rec_fnct_t FaceValuesReconstructFunctor(base_t::m_recEn, m_stencilVals,
    // 							     uMinusHalfNeg, uMinusHalfPos,
    // 							     uPlusHalfNeg,  uPlusHalfPos);

    // EigenJacobianThreeDofPerCellFunctor<> CellJacobianFunctor(jacobian);

    // // auto svals2 = u_stencilVals;
    // // edge_rec_t uMinusHalfNegFD(numDofPerCell);
    // // edge_rec_t uMinusHalfPosFD(numDofPerCell);
    // // edge_rec_t uPlusHalfNegFD(numDofPerCell);
    // // edge_rec_t uPlusHalfPosFD(numDofPerCell);
    // // typename base_t::rec_fnct_t FaceValuesReconstructFunctor2(base_t::m_recEn, svals2,
    // // 					       uMinusHalfNegFD, uMinusHalfPosFD,
    // // 					       uPlusHalfNegFD,  uPlusHalfPosFD);

    // // -------------------
    // // 5. loop
    // // -------------------
    // m_tripletList.clear();
    // const auto dxInv = m_meshObj.dxInv();
    // const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    // const auto & graph = m_meshObj.graph();
    // for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
    //   {
    // 	const auto vIndexDensityDof = smPt*numDofPerCell;

    // 	// the stencil filler populates m_stencilVals with state values
    // 	// for the stencil needed at the current sample mesh point
    // 	FillStencilValuesFunctor(smPt);

    // 	// now that the stencil is ready, do reconstruction which will
    // 	// compute uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
    // 	FaceValuesReconstructFunctor();

    // 	// once reconstruction is done, compute fluxes
    // 	switch(base_t::m_fluxEn)
    // 	{
    // 	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
    // 	  ee::impl::eeRusanovFluxThreeDof(FL, uMinusHalfNeg, uMinusHalfPos, m_gamma);
    // 	  ee::impl::eeRusanovFluxThreeDof(FR, uPlusHalfNeg,  uPlusHalfPos,  m_gamma);

    // 	  ee::impl::eeRusanovFluxJacobianThreeDof(JLneg, JLpos, uMinusHalfNeg, uMinusHalfPos, m_gamma);
    // 	  ee::impl::eeRusanovFluxJacobianThreeDof(JRneg, JRpos, uPlusHalfNeg,  uPlusHalfPos, m_gamma);
    // 	  break;
    // 	}

    // 	CellJacobianFunctor(smPt);

    // 	// store into velocity/rhs
    // 	V(vIndexDensityDof)   = dxInv*(FL(0) - FR(0));
    // 	V(vIndexDensityDof+1) = dxInv*(FL(1) - FR(1));
    // 	V(vIndexDensityDof+2) = dxInv*(FL(2) - FR(2));


    // 	// u_stencilVals = [ u_i-2 (density, mom, e),
    // 	//		     u_i-1 (density, mom, e),
    // 	//		     u_i   (density, mom, e),
    // 	//		     u_i+1 (density, mom, e),
    // 	//		     u_i+2 (density, mom, e)]


    // 	// // the graph gives me the cell IDs wrt to stencil mesh
    // 	// const auto cellGIDStencilMesh      = graph(smPt, 0);
    // 	// const auto leftCellGIDStencilMesh  = graph(smPt, 1);
    // 	// const auto rightCellGIDStencilMesh = graph(smPt, 2);

    // 	// if (base_t::m_probEn == ::pressiodemoapps::Euler1d::PeriodicSmooth)
    // 	// {
    // 	//   addFirstOrderJacobianBlockElementsNoSpecialTreatment(vIndexDensityDof,
    // 	// 						       leftCellGIDStencilMesh*numDofPerCell,
    // 	// 						       cellGIDStencilMesh*numDofPerCell,
    // 	// 						       rightCellGIDStencilMesh*numDofPerCell,
    // 	// 						       JLneg, JLpos, JRneg, JRpos);
    // 	// }
    // 	// else if (base_t::m_probEn == ::pressiodemoapps::Euler1d::Sod or
    // 	// 	 base_t::m_probEn == ::pressiodemoapps::Euler1d::Lax)
    // 	// {
    // 	//   // if we are in inner cells, no special treatment required
    // 	//   if (leftCellGIDStencilMesh != -1 && rightCellGIDStencilMesh != -1){
    // 	//     addFirstOrderJacobianBlockElementsNoSpecialTreatment(vIndexDensityDof,
    // 	// 							 leftCellGIDStencilMesh*numDofPerCell,
    // 	// 							 cellGIDStencilMesh*numDofPerCell,
    // 	// 							 rightCellGIDStencilMesh*numDofPerCell,
    // 	// 							 JLneg, JLpos, JRneg, JRpos);
    // 	//   }

    // 	//   if (leftCellGIDStencilMesh == -1)
    // 	//   {
    // 	//     addFirstOrderJacobianBlockElementsForCellNearLeftBdWithNeumannBC(vIndexDensityDof,
    // 	// 								     cellGIDStencilMesh*numDofPerCell,
    // 	// 								     rightCellGIDStencilMesh*numDofPerCell,
    // 	// 								     JLneg, JLpos, JRneg, JRpos);
    // 	//   }

    // 	//   if (rightCellGIDStencilMesh == -1)
    // 	//   {
    // 	//     addFirstOrderJacobianBlockElementsForCellNearRightBdWithNeumannBC(vIndexDensityDof,
    // 	// 								      rightCellGIDStencilMesh*numDofPerCell,
    // 	// 								      cellGIDStencilMesh*numDofPerCell,
    // 	// 								      JLneg, JLpos, JRneg, JRpos);
    // 	//   }
    //   }

    // m_jacobian.setFromTriplets(m_tripletList.begin(), m_tripletList.end());
    // m_jacobian *= dxInv;
  }

private:
  // if true, Jacobian and velocity are computed together
  bool m_fuseVelocityAndJacobian = false;

  using base_t::m_meshObj;
  using base_t::m_ghostLeft;
  using base_t::m_ghostRight;
  using base_t::m_stencilVals;
  using base_t::m_gamma;
  mutable jacobian_type m_jacobian = {};
  mutable std::vector<Tr> m_tripletList;
};
#endif

}}
#endif
