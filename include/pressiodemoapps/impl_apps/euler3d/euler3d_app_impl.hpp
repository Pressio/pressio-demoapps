
#ifndef PRESSIODEMOAPPS_EULER3D_HPP_
#define PRESSIODEMOAPPS_EULER3D_HPP_

// note that the code below is pretty ugly for now, but it works.
// this is inside impl namespace for a reason, and will need to be
// improved later on but we have a starting point.

#include "../eulerCommon/energy.hpp"
#include "../eulerCommon/fluxes.hpp"
#include "./initial_condition.hpp"
#include "./ghost_filler.hpp"
#include "../stencil_filler.hpp"
#include "../reconstructor_from_stencil.hpp"
#include "../reconstructor_from_state.hpp"

namespace pressiodemoapps{ namespace ee{ namespace impl{

template<
  class scalar_t,
  class mesh_t,
  class state_t,
  class velo_t,
  class ghost_t
  >
class Euler3dAppT
{
  static_assert
  (std::is_same<scalar_t, typename mesh_t::scalar_t>::value, "");

public:
  using index_t		 = typename mesh_t::index_t;
  using scalar_type	 = scalar_t;
  using state_type	 = state_t;
  using velocity_type	 = state_t;
  using ghost_type	 = ghost_t;
  using stencil_values_t = state_type;
  using flux_t		 = state_type;
  using edge_rec_t	 = state_type;

  static constexpr int dimensionality{3};
  static constexpr index_t numDofPerCell{5};
  using normal_vec_t = std::array<scalar_type, dimensionality>;

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
    m_numDofStencilMesh = m_meshObj.stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * numDofPerCell;

    const auto stencilSize = reconstructionTypeToStencilSize(recEn);
    ::pressiodemoapps::resize(m_stencilVals,  numDofPerCell*stencilSize);
    allocateGhosts();
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
    m_numDofStencilMesh = m_meshObj.stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * numDofPerCell;

    const auto stencilSize = reconstructionTypeToStencilSize(recEn);
    ::pressiodemoapps::resize(m_stencilVals,  numDofPerCell*stencilSize);
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
    flux_t FU(numDofPerCell);
    flux_t FD(numDofPerCell);
    flux_t FF(numDofPerCell);
    flux_t FB(numDofPerCell);
    edge_rec_t uMinusHalfNeg(numDofPerCell);
    edge_rec_t uMinusHalfPos(numDofPerCell);
    edge_rec_t uPlusHalfNeg (numDofPerCell);
    edge_rec_t uPlusHalfPos (numDofPerCell);

    fillGhosts(state);
    velocityCellsNearBdImpl(state, timeValue,
			    V, FL, FR, FB, FF, FD, FU,
			    uMinusHalfNeg, uMinusHalfPos,
			    uPlusHalfNeg,  uPlusHalfPos);
    velocityInnerCellsImpl(state, timeValue,
			   V, FL, FR, FB, FF, FD, FU,
			   uMinusHalfNeg, uMinusHalfPos,
			   uPlusHalfNeg,  uPlusHalfPos);
  }

private:
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

    using sfiller_t = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_values_t, state_type, mesh_t, ghost_t>;
    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromStencil<
      edge_rec_t, stencil_values_t>;

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

    const auto dxInv   = m_meshObj.dxInv();
    const auto dyInv   = m_meshObj.dyInv();
    const auto dzInv   = m_meshObj.dzInv();
    auto vEval = [&](index_t vIndex){
      V(vIndex)   = dxInv*(FL(0) - FR(0)) + dyInv*(FB(0) - FF(0)) + dzInv*(FD(0) - FU(0));
      V(vIndex+1) = dxInv*(FL(1) - FR(1)) + dyInv*(FB(1) - FF(1)) + dzInv*(FD(1) - FU(1));
      V(vIndex+2) = dxInv*(FL(2) - FR(2)) + dyInv*(FB(2) - FF(2)) + dzInv*(FD(2) - FU(2));
      V(vIndex+3) = dxInv*(FL(3) - FR(3)) + dyInv*(FB(3) - FF(3)) + dzInv*(FD(3) - FU(3));
      V(vIndex+4) = dxInv*(FL(4) - FR(4)) + dyInv*(FB(4) - FF(4)) + dzInv*(FD(4) - FU(4));
    };

    const auto & graph = m_meshObj.graph();
    const auto & specialRows = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<specialRows.size(); ++it)
    {
      const auto smPt    = specialRows[it];
      const auto cellGID = graph(smPt, 0);

      // X
      StencilFillerX(smPt, it);
      Reconstructor.template operator()<numDofPerCell>();
      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFiveDof(FL, uMinusHalfNeg, uMinusHalfPos, normalX_, m_gamma);
	  eeRusanovFluxFiveDof(FR, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gamma);
	  break;
	}

      // Y
      StencilFillerY(smPt, it);
      Reconstructor.template operator()<numDofPerCell>();
      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFiveDof(FB, uMinusHalfNeg, uMinusHalfPos, normalY_, m_gamma);
	  eeRusanovFluxFiveDof(FF, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gamma);
	  break;
	}

      // Z
      StencilFillerZ(smPt, it);
      Reconstructor.template operator()<numDofPerCell>();
      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFiveDof(FD, uMinusHalfNeg, uMinusHalfPos, normalZ_, m_gamma);
	  eeRusanovFluxFiveDof(FU, uPlusHalfNeg,  uPlusHalfPos,  normalZ_, m_gamma);
	  break;
	}

      const auto vIndex = smPt*numDofPerCell;
      vEval(vIndex);
    }
  }

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

    auto vEval = [&](index_t vIndex){
      V(vIndex)   = dxInv*(FL(0) - FR(0)) + dyInv*(FB(0) - FF(0)) + dzInv*(FD(0) - FU(0));
      V(vIndex+1) = dxInv*(FL(1) - FR(1)) + dyInv*(FB(1) - FF(1)) + dzInv*(FD(1) - FU(1));
      V(vIndex+2) = dxInv*(FL(2) - FR(2)) + dyInv*(FB(2) - FF(2)) + dzInv*(FD(2) - FU(2));
      V(vIndex+3) = dxInv*(FL(3) - FR(3)) + dyInv*(FB(3) - FF(3)) + dzInv*(FD(3) - FU(3));
      V(vIndex+4) = dxInv*(FL(4) - FR(4)) + dyInv*(FB(4) - FF(4)) + dzInv*(FD(4) - FU(4));
    };

    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromState<
      dimensionality, edge_rec_t, state_type, mesh_t>;

    rec_fnct_t ReconstructorX(xAxis, m_recEn, U, m_meshObj,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);

    rec_fnct_t ReconstructorY(yAxis, m_recEn, U, m_meshObj,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);

    rec_fnct_t ReconstructorZ(zAxis, m_recEn, U, m_meshObj,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);

    // deal with cells away from boundaries
    const auto & graph = m_meshObj.graph();
    const auto & rowsIn = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (int it=0; it<rowsIn.size(); ++it)
    {
      const auto smPt = rowsIn[it];
      const auto cellGID = graph(smPt, 0);

      // X
      ReconstructorX.template operator()<numDofPerCell>(smPt);
      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFiveDof(FL, uMinusHalfNeg, uMinusHalfPos, normalX_, m_gamma);
	  eeRusanovFluxFiveDof(FR, uPlusHalfNeg,  uPlusHalfPos,  normalX_, m_gamma);
	  break;
	}

      // Y
      ReconstructorY.template operator()<numDofPerCell>(smPt);
      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFiveDof(FB, uMinusHalfNeg, uMinusHalfPos, normalY_, m_gamma);
	  eeRusanovFluxFiveDof(FF, uPlusHalfNeg,  uPlusHalfPos,  normalY_, m_gamma);
	  break;
	}

      // Z
      ReconstructorZ.template operator()<numDofPerCell>(smPt);
      switch(m_fluxEn)
	{
	case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
	  eeRusanovFluxFiveDof(FD, uMinusHalfNeg, uMinusHalfPos, normalZ_, m_gamma);
	  eeRusanovFluxFiveDof(FU, uPlusHalfNeg,  uPlusHalfPos,  normalZ_, m_gamma);
	  break;
	}

      const auto vIndex = smPt*numDofPerCell;
      vEval(vIndex);
    }

  }

  void fillGhosts(const state_type & U) const
  {
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    if (m_probEn == ::pressiodemoapps::Euler3d::SedovSymmetry)
    {
      using ghost_filler_t  = ::pressiodemoapps::impl::Ghost3dSedov<state_type, mesh_t, ghost_t>;
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

  const normal_vec_t normalX_{1, 0, 0};
  const normal_vec_t normalY_{0, 1, 0};
  const normal_vec_t normalZ_{0, 0, 1};

};

template<class scalar_t, class mesh_t, class state_t, class velo_t, class ghost_t>
constexpr int Euler3dAppT<scalar_t, mesh_t, state_t, velo_t, ghost_t>::dimensionality;

template<class scalar_t, class mesh_t, class state_t, class velo_t, class ghost_t>
constexpr typename mesh_t::index_t Euler3dAppT<scalar_t, mesh_t, state_t, velo_t, ghost_t>::numDofPerCell;


}}}//end namespace
#endif
