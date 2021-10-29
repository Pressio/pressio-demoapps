
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

  // stencil filler functor type
  using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
    dimensionality, numDofPerCell, stencil_vals_container_type,
    state_type, MeshType, ghost_vals_container_type>;

  // reconstructor functor type
  using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromStencilThreeDofPerCell<
    edge_rec_t, stencil_vals_container_type>;

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
		const scalar_type timeValue,
		velocity_type & V) const
  {
    velocityOnlyImpl(state, timeValue, V);
  }

private:
  void velocityOnlyImpl(const state_type & U,
			const scalar_type timeValue,
			velocity_type & V) const
  {

    flux_t FL(numDofPerCell);
    flux_t FR(numDofPerCell);
    edge_rec_t uMinusHalfNeg(numDofPerCell);
    edge_rec_t uMinusHalfPos(numDofPerCell);
    edge_rec_t uPlusHalfNeg (numDofPerCell);
    edge_rec_t uPlusHalfPos (numDofPerCell);
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

    using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell,
      stencil_vals_container_type, state_type, MeshType, ghost_vals_container_type>;
    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromStencilThreeDofPerCell<
      edge_rec_t, stencil_vals_container_type>;

    sfiller_t StencilFiller(stencilSizeNeeded, U, m_meshObj,
			    m_ghostLeft, m_ghostRight, m_stencilVals);
    rec_fnct_t Reconstructor(m_recEn, m_stencilVals,
			     uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg,  uPlusHalfPos);

    const auto dxInv = m_meshObj.dxInv();
    const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
      {
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
		const scalar_type timeValue,
		velocity_type & V) const
  {
    if (base_t::m_recEn == ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder){
      velocityAndJacobianFirstOrderImpl(state, timeValue, V);
    }
    else{
      throw std::runtime_error("Euler1dApp: Jacobian only supported for first order");
    }
  }

  void jacobian(const state_type & state,
		const scalar_type timeValue,
		jacobian_type & J) const
  {
    // assuming jacobian has been computed above
    J = m_jacobian;
  }

private:

  template<class RowIndexType, class FaceJacType>
  void addFirstOrderJacobianBlockElementsNoSpecialTreatment(RowIndexType rowIndex,
							    index_t c_im1,
							    index_t c_i,
							    index_t c_ip1,
							    FaceJacType JLneg,
							    FaceJacType JLpos,
							    FaceJacType JRneg,
							    FaceJacType JRpos) const
  {
    m_tripletList.push_back( Tr(rowIndex,  c_im1,   JLneg(0,0)) );
    m_tripletList.push_back( Tr(rowIndex,  c_im1+1, JLneg(0,1)) );
    m_tripletList.push_back( Tr(rowIndex,  c_im1+2, JLneg(0,2)) );
    m_tripletList.push_back( Tr(rowIndex,  c_i,     JLpos(0,0)-JRneg(0,0)) );
    m_tripletList.push_back( Tr(rowIndex,  c_i+1,   JLpos(0,1)-JRneg(0,0)) );
    m_tripletList.push_back( Tr(rowIndex,  c_i+2,   JLpos(0,2)-JRneg(0,0)) );
    m_tripletList.push_back( Tr(rowIndex,  c_ip1,   -JRpos(0,0)) );
    m_tripletList.push_back( Tr(rowIndex,  c_ip1+1, -JRpos(0,1)) );
    m_tripletList.push_back( Tr(rowIndex,  c_ip1+2, -JRpos(0,2)) );

    m_tripletList.push_back( Tr(rowIndex+1, c_im1,   JLneg(1,0)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_im1+1, JLneg(1,1)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_im1+2, JLneg(1,2)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_i,   JLpos(1,0)-JRneg(1,0)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_i+1, JLpos(1,1)-JRneg(1,0)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_i+2, JLpos(1,2)-JRneg(1,0)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_ip1,   -JRpos(1,0)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_ip1+1, -JRpos(1,1)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_ip1+2, -JRpos(1,2)) );

    m_tripletList.push_back( Tr(rowIndex+2, c_im1,   JLneg(2,0)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_im1+1, JLneg(2,1)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_im1+2, JLneg(2,2)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_i,   JLpos(2,0)-JRneg(2,0)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_i+1, JLpos(2,1)-JRneg(2,0)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_i+2, JLpos(2,2)-JRneg(2,0)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_ip1,   -JRpos(2,0)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_ip1+1, -JRpos(2,1)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_ip1+2, -JRpos(2,2)) );
  }

  template<class RowIndexType, class FaceJacType>
  void addFirstOrderJacobianBlockElementsForCellNearLeftBdWithNeumannBC(RowIndexType rowIndex,
									index_t c_i,
									index_t c_ip1,
									FaceJacType JLneg,
									FaceJacType JLpos,
									FaceJacType JRneg,
									FaceJacType JRpos) const
  {

    m_tripletList.push_back( Tr(rowIndex,  c_i,     JLneg(0,0) + JLpos(0,0) - JRneg(0,0)) );
    m_tripletList.push_back( Tr(rowIndex,  c_i+1,   JLneg(0,1) + JLpos(0,1) - JRneg(0,0)) );
    m_tripletList.push_back( Tr(rowIndex,  c_i+2,   JLneg(0,2) + JLpos(0,2) - JRneg(0,0)) );
    m_tripletList.push_back( Tr(rowIndex,  c_ip1,   -JRpos(0,0)) );
    m_tripletList.push_back( Tr(rowIndex,  c_ip1+1, -JRpos(0,1)) );
    m_tripletList.push_back( Tr(rowIndex,  c_ip1+2, -JRpos(0,2)) );

    m_tripletList.push_back( Tr(rowIndex+1, c_i,   JLneg(1,0) + JLpos(1,0)-JRneg(1,0)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_i+1, JLneg(1,1) + JLpos(1,1)-JRneg(1,0)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_i+2, JLneg(1,2) + JLpos(1,2)-JRneg(1,0)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_ip1,   -JRpos(1,0)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_ip1+1, -JRpos(1,1)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_ip1+2, -JRpos(1,2)) );

    m_tripletList.push_back( Tr(rowIndex+2, c_i,   JLneg(2,0) + JLpos(2,0)-JRneg(2,0)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_i+1, JLneg(2,1) + JLpos(2,1)-JRneg(2,0)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_i+2, JLneg(2,2) + JLpos(2,2)-JRneg(2,0)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_ip1,   -JRpos(2,0)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_ip1+1, -JRpos(2,1)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_ip1+2, -JRpos(2,2)) );
  }

  template<class RowIndexType, class FaceJacType>
  void addFirstOrderJacobianBlockElementsForCellNearRightBdWithNeumannBC(RowIndexType rowIndex,
									index_t c_im1,
									index_t c_i,
									FaceJacType JLneg,
									FaceJacType JLpos,
									FaceJacType JRneg,
									FaceJacType JRpos) const
  {

    m_tripletList.push_back( Tr(rowIndex,  c_im1,   JLneg(0,0)) );
    m_tripletList.push_back( Tr(rowIndex,  c_im1+1, JLneg(0,1)) );
    m_tripletList.push_back( Tr(rowIndex,  c_im1+2, JLneg(0,2)) );
    m_tripletList.push_back( Tr(rowIndex,  c_i,     -JRpos(0,0) + JLpos(0,0)-JRneg(0,0)) );
    m_tripletList.push_back( Tr(rowIndex,  c_i+1,   -JRpos(0,1) + JLpos(0,1)-JRneg(0,0)) );
    m_tripletList.push_back( Tr(rowIndex,  c_i+2,   -JRpos(0,2) + JLpos(0,2)-JRneg(0,0)) );

    m_tripletList.push_back( Tr(rowIndex+1, c_im1,   JLneg(1,0)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_im1+1, JLneg(1,1)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_im1+2, JLneg(1,2)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_i,   -JRpos(1,0) + JLpos(1,0)-JRneg(1,0)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_i+1, -JRpos(1,1) + JLpos(1,1)-JRneg(1,0)) );
    m_tripletList.push_back( Tr(rowIndex+1, c_i+2, -JRpos(1,2) + JLpos(1,2)-JRneg(1,0)) );

    m_tripletList.push_back( Tr(rowIndex+2, c_im1,   JLneg(2,0)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_im1+1, JLneg(2,1)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_im1+2, JLneg(2,2)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_i,   -JRpos(2,0) + JLpos(2,0)-JRneg(2,0)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_i+1, -JRpos(2,1) + JLpos(2,1)-JRneg(2,0)) );
    m_tripletList.push_back( Tr(rowIndex+2, c_i+2, -JRpos(2,2) + JLpos(2,2)-JRneg(2,0)) );
  }

  void fillGhostsIfNeeded(const state_type & U) const
  {
    const auto stencilSizeNeeded = reconstructionTypeToStencilSize(base_t::m_recEn);

    if (base_t::m_probEn == ::pressiodemoapps::Euler1d::Sod or
	base_t::m_probEn == ::pressiodemoapps::Euler1d::Lax)
    {
      using ghost_filler_t  = ::pressiodemoapps::impl::Ghost1dNeumannFiller<
	numDofPerCell, state_type, MeshType, ghost_vals_container_type>;
      ghost_filler_t ghF(stencilSizeNeeded, U, m_meshObj, m_ghostLeft, m_ghostRight);
      ghF();
    }
  }

  void velocityAndJacobianFirstOrderImpl(const state_type & U,
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

    // // -----------------
    // // 2. fill ghosts
    // // -----------------
    // fillGhostsIfNeeded(U);

    // // -------------------
    // // 4. create functors
    // // -------------------
    // typename base_t::sfiller_t FillStencilValuesFunctor(stencilSizeNeeded, U, m_meshObj,
    // 							m_ghostLeft, m_ghostRight,
    // 							m_stencilVals);
    // typename base_t::rec_fnct_t FaceValuesReconstructFunctor(base_t::m_recEn, m_stencilVals,
    // 							     uMinusHalfNeg, uMinusHalfPos,
    // 							     uPlusHalfNeg,  uPlusHalfPos);

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
    // 	// the stencil filler populates m_stencilVals with state values
    // 	// for the stencil needed at the current sample mesh point
    // 	FillStencilValuesFunctor(smPt);

    // 	// now that the stencil is ready, do reconstruction which will
    // 	// compute uMinusHalfNeg, uMinusHalfPos, uPlusHalfNeg,  uPlusHalfPos
    // 	FaceValuesReconstructFunctor.template operator()<numDofPerCell>();

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

    // 	// u_stencilVals = [ u_i-2 (density, mom, e),
    // 	//		     u_i-1 (density, mom, e),
    // 	//		     u_i   (density, mom, e),
    // 	//		     u_i+1 (density, mom, e),
    // 	//		     u_i+2 (density, mom, e)]

    // 	Matrix<double,-1,-1> q(3, u_stencilVals.size());
    // 	// q(:, 0) is d_uLneg/d_ui-2
    // 	// q(:, 1) is d_uLneg/d_ui-1
    // 	// q(:, 2) is d_uLneg/d_ui
    // 	// q(:, 3) is d_uLpos/d_ui-1
    // 	// q(:, 4) is d_uLpos/d_ui
    // 	// q(:, 5) is d_uLpos/d_ui+1

    // 	// q(:, 6)  is d_uRneg/d_ui-1
    // 	// q(:, 7)  is d_uRneg/d_ui
    // 	// q(:, 8)  is d_uRneg/d_ui+1
    // 	// q(:, 9)  is d_uRpos/d_ui
    // 	// q(:, 10) is d_uRpos/d_ui+1
    // 	// q(:, 11) is d_uRpos/d_ui+2

    // 	const auto vIndexDensityDof = smPt*numDofPerCell;
    // 	const auto w0 = graph(smPt, 1);
    // 	const auto e0 = graph(smPt, 2);
    // 	const auto w1 = graph(smPt, 3);
    // 	const auto e1 = graph(smPt, 4);

    // 	const auto num_pertur = u_stencilVals.size()/numDofPerCell;
    // 	for (int j=0; j<num_perturb; ++j)
    // 	{
    // 	  svals2 = m_stencilVals;
    // 	  // if j=0, perturb first numDofPerCell of stencilVals2
    // 	  // if j=1, perturb second numDofPerCell of stencilVals2
    // 	  // ...

    // 	  // assuming that we have svals2 that contains the perturb values
    // 	  // we need to run the reconstruction
    // 	  Reconstructor2.template operator()<numDofPerCell>();
    // 	  // here, it means that, these are overwritten:
    // 	  // uMinusHalfNegFD
    // 	  // uMinusHalfPosFD
    // 	  // uPlusHalfNegFD
    // 	  // uPlusHalfPosFD

    // 	  if (j==0){
    // 	    q(0,0) = (uMinusHalfNegFD(0)-uMinusHalfNeg(0))/eps;
    // 	    q(1,0) = (uMinusHalfNegFD(1)-uMinusHalfNeg(1))/eps;
    // 	    q(2,0) = (uMinusHalfNegFD(2)-uMinusHalfNeg(2))/eps;

    // 	    auto M = JLneg * q.col(0).asDiagonalMatrix();
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof,   w1*numDofPerCell,   M(0,0)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof,   w1*numDofPerCell+1, M(0,1)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof,   w1*numDofPerCell+2, M(0,2)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof+1, w1*numDofPerCell,   M(1,0)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof+1, w1*numDofPerCell+1, M(1,1)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof+1, w1*numDofPerCell+2, M(1,2)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof+2, w1*numDofPerCell,   M(2,0)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof+2, w1*numDofPerCell+1, M(2,1)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof+2, w1*numDofPerCell+2, M(2,2)) );
    // 	  }

    // 	  if (j==1){
    // 	    q(0,1) = (uMinusHalfNegFD(0)-uMinusHalfNeg(0))/eps;
    // 	    q(1,1) = (uMinusHalfNegFD(1)-uMinusHalfNeg(1))/eps;
    // 	    q(2,1) = (uMinusHalfNegFD(2)-uMinusHalfNeg(2))/eps;
    // 	    auto M1 = JLneg * q.col(1).asDiagonalMatrix();

    // 	    q(0,3) = (uMinusHalfPosFD(0)-uMinusHalfPos(0))/eps;
    // 	    q(1,3) = (uMinusHalfPosFD(1)-uMinusHalfPos(1))/eps;
    // 	    q(2,3) = (uMinusHalfPosFD(2)-uMinusHalfPos(2))/eps;
    // 	    auto M2 = JLpos * q.col(3).asDiagonalMatrix();

    // 	    q(0,6) = (uPlusHalfNegFD(0)-uPlusHalfNeg(0))/eps;
    // 	    q(1,6) = (uPlusHalfNegFD(1)-uPlusHalfNeg(1))/eps;
    // 	    q(2,6) = (uPlusHalfNegFD(2)-uPlusHalfNeg(2))/eps;
    // 	    auto M3 = JRneg * q.col(6).asDiagonalMatrix();

    // 	    auto M = M1+M2+M3;
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof,   w0*numDofPerCell,   M(0,0)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof,   w0*numDofPerCell+1, M(0,1)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof,   w0*numDofPerCell+2, M(0,2)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof+1, w0*numDofPerCell,   M(1,0)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof+1, w0*numDofPerCell+1, M(1,1)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof+1, w0*numDofPerCell+2, M(1,2)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof+2, w0*numDofPerCell,   M(2,0)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof+2, w0*numDofPerCell+1, M(2,1)) );
    // 	    m_tripletList.push_back( Tr(vIndexDensityDof+2, w0*numDofPerCell+2, M(2,2)) );
    // 	  }

    // 	  if (j==2){
    // 	    q(0,2) = (uMinusHalfNegFD(0)-uMinusHalfNeg(0))/eps;
    // 	    q(1,2) = (uMinusHalfNegFD(1)-uMinusHalfNeg(1))/eps;
    // 	    q(2,2) = (uMinusHalfNegFD(2)-uMinusHalfNeg(2))/eps;

    // 	    q(0,4) = (uMinusHalfPosFD(0)-uMinusHalfPos(0))/eps;
    // 	    q(1,4) = (uMinusHalfPosFD(1)-uMinusHalfPos(1))/eps;
    // 	    q(2,4) = (uMinusHalfPosFD(2)-uMinusHalfPos(2))/eps;

    // 	    q(0,7) = (uPlusHalfNegFD(0)-uPlusHalfNeg(0))/eps;
    // 	    q(1,7) = (uPlusHalfNegFD(1)-uPlusHalfNeg(1))/eps;
    // 	    q(2,7) = (uPlusHalfNegFD(2)-uPlusHalfNeg(2))/eps;

    // 	    q(0,9) = (uPlusHalfPosFD(0)-uPlusHalfPos(0))/eps;
    // 	    q(1,9) = (uPlusHalfPosFD(1)-uPlusHalfPos(1))/eps;
    // 	    q(2,9) = (uPlusHalfPosFD(2)-uPlusHalfPos(2))/eps;
    // 	  }

    // 	  if (j==3){
    // 	    q(0,5) = (uMinusHalfPosFD(0)-uMinusHalfPos(0))/eps;
    // 	    q(1,5) = (uMinusHalfPosFD(1)-uMinusHalfPos(1))/eps;
    // 	    q(2,5) = (uMinusHalfPosFD(2)-uMinusHalfPos(2))/eps;

    // 	    q(0,8) = (uPlusHalfNegFD(0)-uPlusHalfNeg(0))/eps;
    // 	    q(1,8) = (uPlusHalfNegFD(1)-uPlusHalfNeg(1))/eps;
    // 	    q(2,8) = (uPlusHalfNegFD(2)-uPlusHalfNeg(2))/eps;

    // 	    q(0,10) = (uPlusHalfPosFD(0)-uPlusHalfPos(0))/eps;
    // 	    q(1,10) = (uPlusHalfPosFD(1)-uPlusHalfPos(1))/eps;
    // 	    q(2,10) = (uPlusHalfPosFD(2)-uPlusHalfPos(2))/eps;
    // 	  }

    // 	  if (j==4){
    // 	    q(0,11) = (uPlusHalfPosFD(0)-uPlusHalfPos(0))/eps;
    // 	    q(1,11) = (uPlusHalfPosFD(1)-uPlusHalfPos(1))/eps;
    // 	    q(2,11) = (uPlusHalfPosFD(2)-uPlusHalfPos(2))/eps;
    // 	  }
    // 	}

    // 	// store into velocity/rhs
    // 	const auto vIndexDensityDof = smPt*numDofPerCell;
    // 	V(vIndexDensityDof)   = dxInv*(FL(0) - FR(0));
    // 	V(vIndexDensityDof+1) = dxInv*(FL(1) - FR(1));
    // 	V(vIndexDensityDof+2) = dxInv*(FL(2) - FR(2));

    // 	// the graph gives me the cell IDs wrt to stencil mesh
    // 	const auto cellGIDStencilMesh      = graph(smPt, 0);
    // 	const auto leftCellGIDStencilMesh  = graph(smPt, 1);
    // 	const auto rightCellGIDStencilMesh = graph(smPt, 2);

    // 	if (base_t::m_probEn == ::pressiodemoapps::Euler1d::PeriodicSmooth)
    // 	{
    // 	  addFirstOrderJacobianBlockElementsNoSpecialTreatment(vIndexDensityDof,
    // 							       leftCellGIDStencilMesh*numDofPerCell,
    // 							       cellGIDStencilMesh*numDofPerCell,
    // 							       rightCellGIDStencilMesh*numDofPerCell,
    // 							       JLneg, JLpos, JRneg, JRpos);
    // 	}
    // 	else if (base_t::m_probEn == ::pressiodemoapps::Euler1d::Sod or
    // 		 base_t::m_probEn == ::pressiodemoapps::Euler1d::Lax)
    // 	{
    // 	  // if we are in inner cells, no special treatment required
    // 	  if (leftCellGIDStencilMesh != -1 && rightCellGIDStencilMesh != -1){
    // 	    addFirstOrderJacobianBlockElementsNoSpecialTreatment(vIndexDensityDof,
    // 								 leftCellGIDStencilMesh*numDofPerCell,
    // 								 cellGIDStencilMesh*numDofPerCell,
    // 								 rightCellGIDStencilMesh*numDofPerCell,
    // 								 JLneg, JLpos, JRneg, JRpos);
    // 	  }

    // 	  if (leftCellGIDStencilMesh == -1)
    // 	  {
    // 	    addFirstOrderJacobianBlockElementsForCellNearLeftBdWithNeumannBC(vIndexDensityDof,
    // 									     cellGIDStencilMesh*numDofPerCell,
    // 									     rightCellGIDStencilMesh*numDofPerCell,
    // 									     JLneg, JLpos, JRneg, JRpos);
    // 	  }

    // 	  if (rightCellGIDStencilMesh == -1)
    // 	  {
    // 	    addFirstOrderJacobianBlockElementsForCellNearRightBdWithNeumannBC(vIndexDensityDof,
    // 									      rightCellGIDStencilMesh*numDofPerCell,
    // 									      cellGIDStencilMesh*numDofPerCell,
    // 									      JLneg, JLpos, JRneg, JRpos);
    // 	  }

    // 	}
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
