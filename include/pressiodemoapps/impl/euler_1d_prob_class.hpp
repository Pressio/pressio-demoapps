
#ifndef PRESSIODEMOAPPS_EULER1D_APP_HPP_
#define PRESSIODEMOAPPS_EULER1D_APP_HPP_

// note that the code below is pretty ugly for now, but it works.
// this is inside impl namespace for a reason, and will need to be
// improved later on but we have a starting point.

#include "euler_fluxes.hpp"
#include "euler_flux_jacobian.hpp"
#include "euler_rankine_hugoniot.hpp"
#include "euler_1d_initial_condition.hpp"
#include "euler_1d_ghost_filler.hpp"
#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"
#include "euler_cell_jacobian_first_order.hpp"

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
#include "Eigen/Sparse"
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#include <omp.h>
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

    fillGhostsIfNeeded(state, currentTime);
    velocityCellsNearBdImpl(state, currentTime, V);
    velocityInnerCellsImpl(state, currentTime, V);
  }

protected:
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

private:
  void velocityInnerCellsImpl(const state_type & U,
			      const scalar_type currentTime,
			      velocity_type & V) const
  {
    /*
      for inner cells, we do not need to worry about boundaries,
      we can use the state directly to do the reconstruction.
      we do not need to fill a temporary stencil array values.
    */

    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromState<
      dimensionality, edge_rec_t, state_type, MeshType>;

    flux_t FL(numDofPerCell);
    flux_t FR(numDofPerCell);
    edge_rec_t uMinusHalfNeg(numDofPerCell);
    edge_rec_t uMinusHalfPos(numDofPerCell);
    edge_rec_t uPlusHalfNeg (numDofPerCell);
    edge_rec_t uPlusHalfPos (numDofPerCell);

    rec_fnct_t Reconstructor(m_recEn, U, m_meshObj,
			     uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg,  uPlusHalfPos);

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
			       velocity_type & V) const
  {
    /*
      for cells that are near boundaries, we need use the stencil
      values to handle properly the ghost cells.
      We fill the stencil values needed, then do reconstruction and fluxes.
     */

    flux_t FL(numDofPerCell);
    flux_t FR(numDofPerCell);
    edge_rec_t uMinusHalfNeg(numDofPerCell);
    edge_rec_t uMinusHalfPos(numDofPerCell);
    edge_rec_t uPlusHalfNeg (numDofPerCell);
    edge_rec_t uPlusHalfPos (numDofPerCell);

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
class EigenEuler1dAppWithJacobian
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

public:
  template<class ...Args>
  EigenEuler1dAppWithJacobian(Args && ... args)
    : base_t(std::forward<Args>(args)...),
      m_jacobian(m_numDofSampleMesh, m_numDofStencilMesh)
  {
    initializeJacobian();
  }

  void velocity(const state_type & state,
		const scalar_type currentTime,
		velocity_type & V) const
  {
    base_t::fillGhostsIfNeeded(state, currentTime);

    // jac at left face wrt to left,right reconstructed states
    flux_jac_type JLneg, JLpos;
    // jac at right face wrt to left,right reconstructed states
    flux_jac_type JRneg, JRpos;

    flux_t FL(numDofPerCell);
    flux_t FR(numDofPerCell);
    edge_rec_t uMinusHalfNeg(numDofPerCell);
    edge_rec_t uMinusHalfPos(numDofPerCell);
    edge_rec_t uPlusHalfNeg (numDofPerCell);
    edge_rec_t uPlusHalfPos (numDofPerCell);

    velocityAndJacInnerCellsImpl(state, currentTime,
				 V, FL, FR, JLneg, JLpos, JRneg, JRpos,
				 uMinusHalfNeg, uMinusHalfPos,
				 uPlusHalfNeg,  uPlusHalfPos);

    velocityAndJacNearBDCellsImpl(state, currentTime,
				  V, FL, FR, JLneg, JLpos, JRneg, JRpos,
				  uMinusHalfNeg, uMinusHalfPos,
				  uPlusHalfNeg,  uPlusHalfPos);
  }

  jacobian_type createJacobian() const {
    return m_jacobian;
  }

  void jacobian(const state_type & state,
		const scalar_type timeValue,
		jacobian_type & J) const
  {
    // assuming jacobian has been computed above
    J = m_jacobian;
  }

private:
  void initializeJacobian()
  {
    initializeJacobianFirstOrder();
    if (!m_jacobian.isCompressed()){
      m_jacobian.makeCompressed();
    }
  }

  void initializeJacobianFirstOrder()
  {
    using Tr = Eigen::Triplet<ScalarType>;
    std::vector<Tr> trList;

    //NEED to use first touch here!

    const auto & graph = m_meshObj.graph();
    for (int cell=0; cell<m_meshObj.sampleMeshSize(); ++cell)
      {
	const auto jacRowOfCurrentCellDensity = cell*numDofPerCell;
	const auto ci0  = graph(cell, 0)*numDofPerCell;

	trList.push_back( Tr(jacRowOfCurrentCellDensity,   ci0,   0.0) );
	trList.push_back( Tr(jacRowOfCurrentCellDensity+1, ci0+1, 0.0) );
	trList.push_back( Tr(jacRowOfCurrentCellDensity+2, ci0+2, 0.0) );

	const auto L0 = graph(cell, 1);
	const auto R0 = graph(cell, 2);
	if (L0 != -1){
	  const auto ciL0 = L0*numDofPerCell;
	  trList.push_back( Tr(jacRowOfCurrentCellDensity,   ciL0,   0.0) );
	  trList.push_back( Tr(jacRowOfCurrentCellDensity+1, ciL0+1, 0.0) );
	  trList.push_back( Tr(jacRowOfCurrentCellDensity+2, ciL0+2, 0.0) );
	}

	if (R0 != -1){
	  const auto ciR0 = R0*numDofPerCell;
	  trList.push_back( Tr(jacRowOfCurrentCellDensity,   ciR0,   0.0) );
	  trList.push_back( Tr(jacRowOfCurrentCellDensity+1, ciR0+1, 0.0) );
	  trList.push_back( Tr(jacRowOfCurrentCellDensity+2, ciR0+2, 0.0) );
	}
      }

    m_jacobian.setFromTriplets(trList.begin(), trList.end());
  }

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

    // reconstruct functor
    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromState<
      base_t::dimensionality, edge_rec_t, state_type, MeshType>;
    rec_fnct_t Reconstructor(m_recEn, U, m_meshObj,
			     uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg,  uPlusHalfPos);

    // cell jacobian functor
    using jac_fnct_t = ee::impl::FirstOrderInnerCellJacobianFunctor<
      base_t::dimensionality, numDofPerCell,
      scalar_type, jacobian_type, flux_jac_type, MeshType>;
    jac_fnct_t CellJacobianFunctor(m_jacobian, m_meshObj, JLneg, JLpos, JRneg, JRpos);

    // loop
    const auto dxInv = m_meshObj.dxInv();
    const auto & targetGraphRows = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (int it=0; it<targetGraphRows.size(); ++it)
      {
	const auto smPt = targetGraphRows[it];
	Reconstructor.template operator()<numDofPerCell>(smPt);

	switch(m_fluxEn)
	  {
	  case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	    ee::impl::eeRusanovFluxThreeDof(FL, uMinusHalfNeg, uMinusHalfPos, m_gamma);
	    ee::impl::eeRusanovFluxThreeDof(FR, uPlusHalfNeg,  uPlusHalfPos,  m_gamma);

	    ee::impl::eeRusanovFluxJacobianThreeDof(JLneg, JLpos,
						    uMinusHalfNeg, uMinusHalfPos, m_gamma);
	    ee::impl::eeRusanovFluxJacobianThreeDof(JRneg, JRpos,
						    uPlusHalfNeg, uPlusHalfPos, m_gamma);
	    break;
	  }

	const auto vIndex = smPt*numDofPerCell;
	V(vIndex)   = dxInv*(FL(0) - FR(0));
	V(vIndex+1) = dxInv*(FL(1) - FR(1));
	V(vIndex+2) = dxInv*(FL(2) - FR(2));

	CellJacobianFunctor(smPt);
      }
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

    // stencil filler functor
    const auto stencilSizeNeeded = reconstructionTypeToStencilSize(base_t::m_recEn);
    using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
      base_t::dimensionality, numDofPerCell, stencil_vals_container_type,
      state_type, MeshType, ghost_vals_container_type>;
    sfiller_t FillStencilValuesFunctor(stencilSizeNeeded, U, m_meshObj,
				       m_ghostLeft, m_ghostRight,
				       m_stencilVals);

    // reconstructor functor
    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromStencilThreeDofPerCell<
      edge_rec_t, stencil_vals_container_type>;
    rec_fnct_t FaceValuesReconstructFunctor(base_t::m_recEn, m_stencilVals,
					    uMinusHalfNeg, uMinusHalfPos,
					    uPlusHalfNeg,  uPlusHalfPos);

    // cell jacobian functor
    using jac_fnct_t = ee::impl::FirstOrderNearBDCellJacobianFunctor<
      base_t::dimensionality, numDofPerCell, scalar_type, jacobian_type, flux_jac_type, MeshType>;
    jac_fnct_t CellJacobianFunctor(m_jacobian, m_meshObj, JLneg, JLpos, JRneg, JRpos);

    // loop over cells
    const auto dxInv = m_meshObj.dxInv();
    const auto & targetGraphRows = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<targetGraphRows.size(); ++it)
    {
      const auto smPt = targetGraphRows[it];
      const auto vIndexCurrentCellDensity = smPt*numDofPerCell;

      FillStencilValuesFunctor(smPt);
      FaceValuesReconstructFunctor();

      switch(m_fluxEn)
      {
      case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	ee::impl::eeRusanovFluxThreeDof(FL, uMinusHalfNeg, uMinusHalfPos, m_gamma);
	ee::impl::eeRusanovFluxThreeDof(FR, uPlusHalfNeg,  uPlusHalfPos,  m_gamma);

	ee::impl::eeRusanovFluxJacobianThreeDof(JLneg, JLpos, uMinusHalfNeg, uMinusHalfPos, m_gamma);
	ee::impl::eeRusanovFluxJacobianThreeDof(JRneg, JRpos, uPlusHalfNeg,  uPlusHalfPos, m_gamma);
	break;
      }

      V(vIndexCurrentCellDensity)   = dxInv*(FL(0) - FR(0));
      V(vIndexCurrentCellDensity+1) = dxInv*(FL(1) - FR(1));
      V(vIndexCurrentCellDensity+2) = dxInv*(FL(2) - FR(2));

      CellJacobianFunctor(smPt);
    }
  }

private:
  using base_t::m_meshObj;
  using base_t::m_ghostLeft;
  using base_t::m_ghostRight;
  using base_t::m_stencilVals;
  using base_t::m_gamma;
  using base_t::m_numDofStencilMesh;
  using base_t::m_numDofSampleMesh;
  using base_t::m_recEn;
  using base_t::m_fluxEn;
  mutable jacobian_type m_jacobian = {};
};
#endif

}}
#endif
