
#ifndef PRESSIODEMOAPPS_EULER2D_HPP_
#define PRESSIODEMOAPPS_EULER2D_HPP_

// note that the code below is pretty ugly for now, but it works.
// this is inside impl namespace for a reason, and will need to be
// improved later on but we have a starting point.

namespace pressiodemoapps{ namespace ee{ namespace impl{

template<
  class scalar_t,
  class mesh_t,
  class state_t,
  class velo_t,
  class ghost_t
  >
class Euler2dAppT
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

  static constexpr int dimensionality{2};
  static constexpr index_t numDofPerCell{4};

public:
#if not defined PRESSIODEMOAPPS_ENABLE_BINDINGS
  Euler2dAppT(const mesh_t & meshObj,
	      pressiodemoapps::reconstructionEnum recEn,
	      pressiodemoapps::euler2dproblemsEnum probEn,
	      int icIdentifier = 1)
    : m_recEn(recEn), m_probEn(probEn), m_icIdentifier(icIdentifier), m_meshObj(meshObj)
  {
    m_numDofStencilMesh = m_meshObj.stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * numDofPerCell;

    const auto stencilSize = reconstructionEnumToStencilSize(recEn);
    ::pressiodemoapps::resize(m_stencilVals,  numDofPerCell*stencilSize);
    allocateGhosts();
  }

#else
  // note that when doing bindings, I need to first construct
  // ghost with {1,1} just so that the numpy array picks up they
  // are 2dim array otherwise it thinks they are 1d array.
  // The right allocation for these is then done inside allocateGhosts.
  Euler2dAppT(const mesh_t & meshObj,
	      pressiodemoapps::reconstructionEnum recEn,
	      pressiodemoapps::euler2dproblemsEnum probEn,
	      int icIdentifier = 1)
    : m_recEn(recEn), m_probEn(probEn), m_icIdentifier(icIdentifier),
      m_meshObj(meshObj),
      m_stencilVals(1),
      m_ghostLeft({1,1}),
      m_ghostTop({1,1}),
      m_ghostRight({1,1}),
      m_ghostBottom({1,1})
  {
    m_numDofStencilMesh = m_meshObj.stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * numDofPerCell;

    const auto stencilSize = reconstructionEnumToStencilSize(recEn);
    ::pressiodemoapps::resize(m_stencilVals,  numDofPerCell*stencilSize);
    allocateGhosts();
  }
#endif

  state_type initialCondition() const
  {
    state_type IC(m_numDofStencilMesh);

    switch(m_probEn)
      {
      case pressiodemoapps::euler2dproblemsEnum::periodic:{
	return IC;
      }

      case pressiodemoapps::euler2dproblemsEnum::sedov:{
	sedov2dInitialCondition(IC, m_meshObj, numDofPerCell, m_gamma);
	return IC;
      }

      case pressiodemoapps::euler2dproblemsEnum::riemann:{
	if( m_icIdentifier == 1){
	  riemann2dInitialCondition1(IC, m_meshObj, numDofPerCell, m_gamma);
	  return IC;
	}
	else if (m_icIdentifier == 2){
	  riemann2dInitialCondition2(IC, m_meshObj, numDofPerCell, m_gamma);
	  return IC;
	}
	else{
	  throw std::runtime_error("invalid IC");
	}
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
    edge_rec_t uMinusHalfNeg(numDofPerCell);
    edge_rec_t uMinusHalfPos(numDofPerCell);
    edge_rec_t uPlusHalfNeg (numDofPerCell);
    edge_rec_t uPlusHalfPos (numDofPerCell);

    fillGhosts(state);
    velocityCellsNearBdImpl(state, timeValue, V, FL, FR, FD, FU,
			    uMinusHalfNeg, uMinusHalfPos,
			    uPlusHalfNeg,  uPlusHalfPos);
    velocityInnerCellsImpl(state, timeValue, V, FL, FR, FD, FU,
			   uMinusHalfNeg, uMinusHalfPos,
			   uPlusHalfNeg,  uPlusHalfPos);
  }

private:
  void velocityCellsNearBdImpl(const state_type & U,
			       const scalar_type timeIn,
			       velocity_type & V,
			       flux_t & FL,
			       flux_t & FR,
			       flux_t & FD,
			       flux_t & FU,
			       edge_rec_t & uMinusHalfNeg,
			       edge_rec_t & uMinusHalfPos,
			       edge_rec_t & uPlusHalfNeg,
			       edge_rec_t & uPlusHalfPos) const
  {

    using sfiller_t = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_values_t, state_type, mesh_t, ghost_t>;
    using rec_fnct_t =
      ::pressiodemoapps::impl::ReconstructorFromStencil<edge_rec_t, stencil_values_t>;

    const auto stencilSize = reconstructionEnumToStencilSize(m_recEn);
    sfiller_t StencilFillerX(stencilSize, U, m_meshObj,
			     m_ghostLeft, m_ghostRight,
			     m_stencilVals, 1);

    sfiller_t StencilFillerY(stencilSize, U, m_meshObj,
			     m_ghostBottom, m_ghostTop,
			     m_stencilVals, 2);

    rec_fnct_t Reconstructor(m_recEn, m_stencilVals,
			     uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg,  uPlusHalfPos);

    const auto dxInv   = m_meshObj.dxInv();
    const auto dyInv   = m_meshObj.dyInv();
    auto vEval = [&](index_t vIndex){
      V(vIndex)   = dxInv*(FL(0) - FR(0)) + dyInv*(FD(0) - FU(0));
      V(vIndex+1) = dxInv*(FL(1) - FR(1)) + dyInv*(FD(1) - FU(1));
      V(vIndex+2) = dxInv*(FL(2) - FR(2)) + dyInv*(FD(2) - FU(2));
      V(vIndex+3) = dxInv*(FL(3) - FR(3)) + dyInv*(FD(3) - FU(3));
    };

    const auto & graph = m_meshObj.graph();
    const auto & specialRows = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<specialRows.size(); ++it)
    {
      const auto smPt    = specialRows[it];
      const auto cellGID = graph(smPt, 0);
      const auto uIndex  = cellGID*numDofPerCell;

      // X
      StencilFillerX(smPt, it);
      Reconstructor.template operator()<numDofPerCell>();
      eeRusanovFluxFourDof(FL, uMinusHalfNeg,
			   uMinusHalfPos, normalX_, m_gamma);
      eeRusanovFluxFourDof(FR, uPlusHalfNeg,
			   uPlusHalfPos,  normalX_, m_gamma);

      // Y
      StencilFillerY(smPt, it);
      Reconstructor.template operator()<numDofPerCell>();
      eeRusanovFluxFourDof(FD, uMinusHalfNeg,
			   uMinusHalfPos, normalY_, m_gamma);
      eeRusanovFluxFourDof(FU, uPlusHalfNeg,
			   uPlusHalfPos,  normalY_, m_gamma);

      const auto vIndex = smPt*numDofPerCell;
      vEval(vIndex);
    }
  }

  void velocityInnerCellsImpl(const state_type & U,
			      const scalar_type timeIn,
			      velocity_type & V,
			      flux_t & FL,
			      flux_t & FR,
			      flux_t & FD,
			      flux_t & FU,
			      edge_rec_t & uMinusHalfNeg,
			      edge_rec_t & uMinusHalfPos,
			      edge_rec_t & uPlusHalfNeg,
			      edge_rec_t & uPlusHalfPos) const
  {

    const auto dxInv   = m_meshObj.dxInv();
    const auto dyInv   = m_meshObj.dyInv();
    auto vEval = [&](index_t vIndex){
      V(vIndex)   = dxInv*(FL(0) - FR(0)) + dyInv*(FD(0) - FU(0));
      V(vIndex+1) = dxInv*(FL(1) - FR(1)) + dyInv*(FD(1) - FU(1));
      V(vIndex+2) = dxInv*(FL(2) - FR(2)) + dyInv*(FD(2) - FU(2));
      V(vIndex+3) = dxInv*(FL(3) - FR(3)) + dyInv*(FD(3) - FU(3));
    };

    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromState<
      edge_rec_t, state_type, mesh_t>;

    rec_fnct_t ReconstructorX(1, m_recEn, U, m_meshObj,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);

    rec_fnct_t ReconstructorY(2, m_recEn, U, m_meshObj,
			      uMinusHalfNeg, uMinusHalfPos,
			      uPlusHalfNeg,  uPlusHalfPos);

    // deal with cells away from boundaries
    const auto & graph = m_meshObj.graph();
    const auto & rowsIn = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (int it=0; it<rowsIn.size(); ++it){
      const auto smPt = rowsIn[it];
      const auto cellGID = graph(smPt, 0);
      const auto uIndex = cellGID*numDofPerCell;

      // X
      ReconstructorX.template operator()<numDofPerCell>(smPt);
      eeRusanovFluxFourDof(FL, uMinusHalfNeg,
			   uMinusHalfPos, normalX_, m_gamma);
      eeRusanovFluxFourDof(FR, uPlusHalfNeg,
			   uPlusHalfPos,  normalX_, m_gamma);

      // Y
      ReconstructorY.template operator()<numDofPerCell>(smPt);
      eeRusanovFluxFourDof(FD, uMinusHalfNeg,
			   uMinusHalfPos, normalY_, m_gamma);
      eeRusanovFluxFourDof(FU, uPlusHalfNeg,
			   uPlusHalfPos,  normalY_, m_gamma);

      const auto vIndex = smPt*numDofPerCell;
      vEval(vIndex);
    }

  }

  void fillGhosts(const state_type & U) const
  {
    const auto stencilSize = reconstructionEnumToStencilSize(m_recEn);
    if (m_probEn == pressiodemoapps::euler2dproblemsEnum::sedov or
	m_probEn == pressiodemoapps::euler2dproblemsEnum::riemann)
    {
      using ghost_filler_t  = ::pressiodemoapps::impl::Ghost2dNeumannFiller<
	numDofPerCell, state_type, mesh_t, ghost_t>;
      ghost_filler_t ghF(stencilSize, U, m_meshObj,
			 m_ghostLeft, m_ghostTop,
			 m_ghostRight, m_ghostBottom);

      const auto & rowsBd = m_meshObj.graphRowsOfCellsNearBd();
      if (stencilSize==3){
	for (int it=0; it<rowsBd.size(); ++it){
	  ghF.template operator()<3>(rowsBd[it], it);
	}
      }else{
	for (int it=0; it<rowsBd.size(); ++it){
	  ghF.template operator()<7>(rowsBd[it], it);
	}
      }
    }
  }

private:
  void allocateGhosts()
  {
    /*
      for stencil = 3, at leftBoundary:
      ---------------
      |	 0,1,2,3   ||
      | rho,       ||
      | rho u,	   ||
      | rho v,     ||
      | E	   ||
      ---------------

      for stencil = 7, at leftBoundary:
      --------------------------------------
      |	 8,9,10,11  | 4,5,6,7 |  0,1,2,3  ||
      |	     	    |
      | rho,        | rho,    | rho       ||
      | rho u,	    | rho*u   | rho*u     ||
      | rho v,      | rho*v   | rho*v     ||
      | E	    | E       | E         ||
      --------------------------------------
     */

    const auto stencilSize    = reconstructionEnumToStencilSize(m_recEn);
    const auto numGhostValues = numDofPerCell*((stencilSize-1)/2);

    const index_t s1 = m_meshObj.numCellsBd();
    pressiodemoapps::resize(m_ghostLeft,   s1, numGhostValues);
    pressiodemoapps::resize(m_ghostTop,    s1, numGhostValues);
    pressiodemoapps::resize(m_ghostRight,  s1, numGhostValues);
    pressiodemoapps::resize(m_ghostBottom, s1, numGhostValues);
  }

private:
  const scalar_type m_gamma{1.4};

  pressiodemoapps::reconstructionEnum m_recEn =
    pressiodemoapps::reconstructionEnum::firstOrder;

  pressiodemoapps::euler2dproblemsEnum m_probEn;

  // which initial condition to use
  int m_icIdentifier = 1;

  const mesh_t & m_meshObj;
  mutable stencil_values_t m_stencilVals;

  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  mutable ghost_t m_ghostLeft;
  mutable ghost_t m_ghostTop;
  mutable ghost_t m_ghostRight;
  mutable ghost_t m_ghostBottom;

  const std::array<scalar_type, 2> normalX_{1, 0};
  const std::array<scalar_type, 2> normalY_{0, 1};

};

}}}//end namespace
#endif
