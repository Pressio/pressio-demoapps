
#ifndef PRESSIODEMOAPPS_DIFF_REACTION_2D_IMPL_HPP_
#define PRESSIODEMOAPPS_DIFF_REACTION_2D_IMPL_HPP_

#include "functor_fill_stencil.hpp"
#include "diffusion_reaction_2d_ghost_filler.hpp"

namespace pressiodemoapps{ namespace impldiffreac{

template<class scalar_t>
struct DefaultSourceF2d
{
  void operator()(const scalar_t & x,
		  const scalar_t & y,
		  const scalar_t & evaltime,
		  scalar_t & value)
  {
    (void) evaltime;
    value = std::sin(M_PI*x*(y-0.2)) * 4.*std::sin(4.*M_PI*y*x);
  }
};

template<
  class scalar_t,
  class mesh_t,
  class state_t,
  class velo_t,
  class ghost_t,
  class jacobian_t = void
  >
class DiffReac2dApp
{

public:
  using index_t		 = typename mesh_t::index_t;
  using stencil_values_t = state_t;
  using scalar_type	 = scalar_t;
  using state_type	 = state_t;
  using velocity_type	 = velo_t;
  using ghost_type	 = ghost_t;

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
  using jacobian_type	 = jacobian_t;
  typedef Eigen::Triplet<scalar_type> Tr;
#endif

  static constexpr int dimensionality{2};
  static constexpr int numDofPerCell{1};

public:

  template<class SourceT>
  DiffReac2dApp(const mesh_t & meshObj,
		::pressiodemoapps::DiffusionReaction2d probEnum,
		::pressiodemoapps::ViscousFluxReconstruction recEnum,
		SourceT sf,
		scalar_t diffusionCoeff,
		scalar_t reactionCoeff)
    : m_meshObj(meshObj), m_probEn(probEnum), m_recEn(recEnum)
  {
    if (m_meshObj.stencilSize() != 3){
      throw std::runtime_error("DiffusionReaction2d currently, only supports 3-pt stencil");
    }

    m_numDofStencilMesh = m_meshObj.stencilMeshSize()*numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize()*numDofPerCell;

    m_sourceFunctor = sf;
    m_diffusionCoeff = diffusionCoeff;
    m_reactionCoeff = reactionCoeff;

    const auto stencilSize = reconstructionTypeToStencilSize(recEnum);
    ::pressiodemoapps::resize(m_stencilVals, numDofPerCell*stencilSize);
    allocateGhosts();

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
    m_jacobian.resize(m_numDofSampleMesh, m_numDofStencilMesh);
#endif
  }

  DiffReac2dApp(const mesh_t & meshObj,
		 ::pressiodemoapps::DiffusionReaction2d probEnum,
		 ::pressiodemoapps::ViscousFluxReconstruction recEnum,
		 scalar_t diffusionCoeff,
		 scalar_t reactionCoeff)
    : DiffReac2dApp(meshObj, probEnum, recEnum,
		     DefaultSourceF2d<scalar_t>(), diffusionCoeff, reactionCoeff)
  {}

  DiffReac2dApp(const mesh_t & meshObj,
		 ::pressiodemoapps::DiffusionReaction2d probEnum,
		 ::pressiodemoapps::ViscousFluxReconstruction recEnum)
    : DiffReac2dApp(meshObj, probEnum, recEnum,
		    DefaultSourceF2d<scalar_t>(), 0.01, 0.01)
  {}

  index_t totalDofSampleMesh()  const{ return m_numDofSampleMesh; }
  index_t totalDofStencilMesh() const{ return m_numDofStencilMesh; }

  state_type initialCondition() const
  {
    state_type ic(m_numDofStencilMesh);
    if (m_probEn == ::pressiodemoapps::DiffusionReaction2d::ProblemA){
      for (int i=0; i<::pressiodemoapps::extent(ic,0); ++i){
	ic(i) = scalar_t(0);
      }
    }
    return ic;
  }

  velocity_type createVelocity() const {
    velocity_type V(m_numDofSampleMesh);
    return V;
  }

  void velocity(const state_type & state,
		const scalar_type time,
		velocity_type & v) const
  {
    fillGhostsIfNeeded(state, time);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
    if (!m_onlyComputeVelocity){
      m_tripletList.clear();
    }
#endif

    velocityCellsNearBdImpl(state, time, v);
    velocityInnerCellsImpl(state, time, v);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
    if (!m_onlyComputeVelocity){
      m_jacobian.setFromTriplets(m_tripletList.begin(), m_tripletList.end());
    }
#endif
  }

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
  jacobian_type createJacobian() const{
    jacobian_type JJ(m_numDofSampleMesh, m_numDofStencilMesh);
    return JJ;
  }

  void jacobian(const state_type & state,
		const scalar_type time,
		jacobian_type & J) const
  {
    if (!m_onlyComputeVelocity){
      // relies on jacobian been computed in velocity
      J = m_jacobian;
    }
  }

  // the Jacobian is by default fused with the velocity,
  // this method allows one to disable the jacobian
  // so only velocity is computed
  void disableJacobian() {
    m_onlyComputeVelocity = true;
  }
#endif

private:
  void fillGhostsIfNeeded(const state_type & U,
			  scalar_type /*currTime*/) const
  {
    if (m_probEn == ::pressiodemoapps::DiffusionReaction2d::ProblemA)
    {
      using ghost_filler_t  = ::pressiodemoapps::impldiffreac::GhostFillerProblemA2d<
	state_type, mesh_t, ghost_t>;

      const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_recEn);
      ghost_filler_t ghF(stencilSizeNeeded, U, m_meshObj,
			 m_ghostLeft, m_ghostFront,
			 m_ghostRight, m_ghostBack);

      const auto & rowsBd = m_meshObj.graphRowsOfCellsNearBd();
      for (int it=0; it<rowsBd.size(); ++it){
	ghF(rowsBd[it], it);
      }
    }
  }

  void velocityCellsNearBdImpl(const state_type & U,
			       const scalar_type t,
			       velocity_type & V) const
  {
    // note that numDofPerCell == 1, so we omit it below

    constexpr int xAxis = 1;
    constexpr int yAxis = 2;

    // stencil filler needed because we are doing cells near boundaries
    using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_values_t, state_type, mesh_t, ghost_t>;

    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    sfiller_t StencilFillerX(stencilSize, U, m_meshObj,
			     m_ghostLeft, m_ghostRight,
			     m_stencilVals, xAxis);
    sfiller_t StencilFillerY(stencilSize, U, m_meshObj,
			     m_ghostBack, m_ghostFront,
			     m_stencilVals, yAxis);

    const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto dyInvSq  = m_meshObj.dyInv()*m_meshObj.dyInv();
    const auto & x      = m_meshObj.viewX();
    const auto & y      = m_meshObj.viewY();
    const auto & graph  = m_meshObj.graph();
    constexpr auto two  = static_cast<scalar_t>(2);
    constexpr auto three= static_cast<scalar_t>(3);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
    const auto twoReacCoeff = m_reactionCoeff*two;
    const auto diffDxInvSq  = m_diffusionCoeff*dxInvSq;
    const auto diffDyInvSq  = m_diffusionCoeff*dyInvSq;
#endif

    const auto & rows   = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<rows.size(); ++it)
    {
      const auto smPt        = rows[it];
      const auto uIndex	     = graph(smPt, 0);
      const auto uIndexLeft  = graph(smPt, 1);
      const auto uIndexFront = graph(smPt, 2);
      const auto uIndexRight = graph(smPt, 3);
      const auto uIndexBack  = graph(smPt, 4);

      // compute source, store into V
      m_sourceFunctor(x(uIndex), y(uIndex), t, V(smPt));

      // add to V reaction contribution
      V(smPt) += m_reactionCoeff*U(uIndex)*U(uIndex);

      // *** add X contribution of diffusion ***
      StencilFillerX(smPt, it);
      V(smPt) += dxInvSq*m_diffusionCoeff*( m_stencilVals(2) - two*m_stencilVals(1) + m_stencilVals(0) );

      // *** add Y contribution of diffusion ***
      StencilFillerY(smPt, it);
      V(smPt) += dyInvSq*m_diffusionCoeff*( m_stencilVals(2) - two*m_stencilVals(1) + m_stencilVals(0) );

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (m_probEn == ::pressiodemoapps::DiffusionReaction2d::ProblemA &&
	  !m_onlyComputeVelocity)
      {
	auto selfValue = -two*diffDxInvSq -two*diffDyInvSq + twoReacCoeff*m_stencilVals(1);

	if (uIndexLeft != -1){

	  m_tripletList.push_back( Tr(smPt, uIndexLeft, diffDxInvSq) );
	}else{
	  selfValue += -diffDxInvSq;
	}

	if (uIndexFront != -1){
	  m_tripletList.push_back( Tr(smPt, uIndexFront, diffDyInvSq) );
	}else{
	  selfValue += -diffDyInvSq;
	}

	if (uIndexRight != -1){
	  m_tripletList.push_back( Tr(smPt, uIndexRight, diffDxInvSq) );
	}
	else{
	  selfValue += -diffDxInvSq;
	}

	if (uIndexBack != -1){
	  m_tripletList.push_back( Tr(smPt, uIndexBack, diffDyInvSq) );
	}
	else{
	  selfValue += -diffDyInvSq;
	}

	// note that this MUST bere here because selfValue takes
	// multiple contributions above
	m_tripletList.push_back( Tr(smPt, uIndex, selfValue) );
      }
#endif
    }
  }

  void velocityInnerCellsImpl(const state_type & U,
			      const scalar_type t,
			      velocity_type & V) const
  {
    // note that numDofPerCell == 1, so we omit it below

    const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto dyInvSq  = m_meshObj.dyInv()*m_meshObj.dyInv();
    const auto & x      = m_meshObj.viewX();
    const auto & y      = m_meshObj.viewY();
    const auto & graph  = m_meshObj.graph();
    constexpr auto two  = static_cast<scalar_t>(2);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
    const auto twoReacCoeff = m_reactionCoeff*two;
    const auto diffDxInvSq  = m_diffusionCoeff*dxInvSq;
    const auto diffDyInvSq  = m_diffusionCoeff*dyInvSq;
#endif

    const auto & rows   = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (std::size_t it=0; it<rows.size(); ++it)
    {
      const auto smPt        = rows[it];
      const auto uIndex	     = graph(smPt, 0);
      const auto uIndexLeft  = graph(smPt, 1);
      const auto uIndexFront = graph(smPt, 2);
      const auto uIndexRight = graph(smPt, 3);
      const auto uIndexBack  = graph(smPt, 4);

      // compute source, store into V
      m_sourceFunctor(x(uIndex), y(uIndex), t, V(smPt));

      // ADD to V reaction contribution
      V(smPt) += m_reactionCoeff*U(uIndex)*U(uIndex);

      // ADD to V x diffusion contribution
      V(smPt) += dxInvSq*m_diffusionCoeff*( U(uIndexRight) - two*U(uIndex) + U(uIndexLeft) );

      // ADD to V x diffusion contribution
      V(smPt) += dyInvSq*m_diffusionCoeff*( U(uIndexFront) - two*U(uIndex) + U(uIndexBack) );

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (m_probEn == ::pressiodemoapps::DiffusionReaction2d::ProblemA &&
	  !m_onlyComputeVelocity)
      {
	const auto jvalueself = -two*diffDxInvSq -two*diffDyInvSq + twoReacCoeff*U(uIndex);
	m_tripletList.push_back( Tr(smPt, uIndex,      jvalueself) );
	m_tripletList.push_back( Tr(smPt, uIndexLeft,  diffDxInvSq) );
	m_tripletList.push_back( Tr(smPt, uIndexFront, diffDyInvSq) );
	m_tripletList.push_back( Tr(smPt, uIndexRight, diffDxInvSq) );
	m_tripletList.push_back( Tr(smPt, uIndexBack,  diffDyInvSq) );
      }
#endif
    }
  }

private:
  void allocateGhosts()
  {
    const auto stencilSize    = reconstructionTypeToStencilSize(m_recEn);
    const auto numGhostValues = numDofPerCell*((stencilSize-1)/2);

    const index_t s1 = m_meshObj.numCellsBd();
    ::pressiodemoapps::resize(m_ghostLeft,  s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostFront, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostBack,  s1, numGhostValues);
  }

private:
  const mesh_t & m_meshObj;
  ::pressiodemoapps::DiffusionReaction2d m_probEn;
  ::pressiodemoapps::ViscousFluxReconstruction m_recEn;

  std::function<void(const scalar_t & /*x*/,
		     const scalar_t & /*y*/,
		     const scalar_t & /*time*/,
		     scalar_t &)> m_sourceFunctor;
  scalar_t m_diffusionCoeff = {};
  scalar_t m_reactionCoeff = {};

  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};
  mutable stencil_values_t m_stencilVals;

  mutable ghost_t m_ghostLeft;
  mutable ghost_t m_ghostFront;
  mutable ghost_t m_ghostRight;
  mutable ghost_t m_ghostBack;

  const std::array<scalar_type, 2> normalX_{1, 0};
  const std::array<scalar_type, 2> normalY_{0, 1};

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
  mutable jacobian_type m_jacobian = {};
  mutable std::vector<Tr> m_tripletList;
  bool m_onlyComputeVelocity = false;
#endif
};

}}
#endif
