
#ifndef PRESSIODEMOAPPS_DIFF_REACTION_1D_IMPL_HPP_
#define PRESSIODEMOAPPS_DIFF_REACTION_1D_IMPL_HPP_

#include "functor_fill_stencil.hpp"
#include "diffusion_reaction_1d_ghost_filler.hpp"

namespace pressiodemoapps{ namespace impldiffreac{

template<class scalar_t>
struct DefaultSourceF
{
  void operator()(const scalar_t & x, const scalar_t & evaltime, scalar_t & value){
    (void) evaltime;
    value = std::sin(M_PI*x) * x*x * 4.*std::cos(4.*M_PI*x);
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
class DiffReac1dApp
{

public:
  using index_t		 = typename mesh_t::index_t;
  using stencil_values_t = state_t;
  using scalar_type	 = scalar_t;
  using state_type	 = state_t;
  using velocity_type	 = state_t;
  using ghost_type	 = ghost_t;

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
  using jacobian_type	 = jacobian_t;
  typedef Eigen::Triplet<scalar_type> Tr;
#endif

  static constexpr int dimensionality{1};
  static constexpr int numDofPerCell{1};

public:

  template<class SourceT>
  DiffReac1dApp(const mesh_t & meshObj,
		::pressiodemoapps::DiffusionReaction1d probEnum,
		::pressiodemoapps::ViscousFluxReconstruction recEnum,
		SourceT sf,
		scalar_t diffusionCoeff,
		scalar_t reactionCoeff)
    : m_meshObj(meshObj), m_probEn(probEnum), m_recEn(recEnum)
  {
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

    // make sure that no matter what, the first and last sample cells
    // are near the boundary. Even if we are doing sample mesh, we need
    // to keep the cells near boundaries.
    const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    const auto & graph = m_meshObj.graph();
    if (graph(0,1) != -1){
      throw std::runtime_error("Point not on boundary, something wrong");
    }
    if (graph(sampleMeshSize-1, 2) != -1){
      throw std::runtime_error("Point not on boundary, something wrong");
    }
  }

  DiffReac1dApp(const mesh_t & meshObj,
		 ::pressiodemoapps::DiffusionReaction1d probEnum,
		 ::pressiodemoapps::ViscousFluxReconstruction recEnum,
		 scalar_t diffusionCoeff,
		 scalar_t reactionCoeff)
    : DiffReac1dApp(meshObj, probEnum, recEnum,
		     DefaultSourceF<scalar_t>(), diffusionCoeff, reactionCoeff)
  {}

  DiffReac1dApp(const mesh_t & meshObj,
		 ::pressiodemoapps::DiffusionReaction1d probEnum,
		 ::pressiodemoapps::ViscousFluxReconstruction recEnum)
    : DiffReac1dApp(meshObj, probEnum, recEnum,
		    DefaultSourceF<scalar_t>(), 0.01, 0.01)
  {}

  index_t totalDofSampleMesh()  const{ return m_numDofSampleMesh; }
  index_t totalDofStencilMesh() const{ return m_numDofStencilMesh; }

  state_type initialCondition() const
  {
    state_type ic(m_numDofStencilMesh);
    if (m_probEn == ::pressiodemoapps::DiffusionReaction1d::ProblemA){
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
    if (m_probEn == ::pressiodemoapps::DiffusionReaction1d::ProblemA)
    {
      using ghost_filler_t  = ::pressiodemoapps::impldiffreac::GhostFillerProblemA1d<
	state_type, mesh_t, ghost_t>;

      const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_recEn);
      ghost_filler_t ghF(stencilSizeNeeded, U, m_meshObj, m_ghostLeft, m_ghostRight);
      ghF();
    }
  }

  void velocityCellsNearBdImpl(const state_type & U,
			       const scalar_type t,
			       velocity_type & V) const
  {
    // note that numDofPerCell == 1, so we omit it below

    // stencil filler needed because we are doing cells near boundaries
    using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_values_t, state_type, mesh_t, ghost_t>;
    sfiller_t StencilFiller(reconstructionTypeToStencilSize(m_recEn),
			    U, m_meshObj, m_ghostLeft, m_ghostRight, m_stencilVals);

    const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto & graph  = m_meshObj.graph();
    const auto & x      = m_meshObj.viewX();
    constexpr auto two  = static_cast<scalar_t>(2);
    constexpr auto three= static_cast<scalar_t>(3);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
    const auto twoReacCoeff = m_reactionCoeff*two;
    const auto diffDxInvSq  = m_diffusionCoeff*dxInvSq;
#endif

    const auto & rows   = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<rows.size(); ++it)
    {
      const auto smPt        = rows[it];
      const auto uIndex	     = graph(smPt, 0);
      const auto uIndexLeft  = graph(smPt, 1);
      const auto uIndexRight = graph(smPt, 2);

      StencilFiller(smPt);

      // compute source, store into V
      m_sourceFunctor(x(uIndex), t, V(smPt));

      // ADD to V reaction contribution
      V(smPt) += m_reactionCoeff*U(uIndex)*U(uIndex);

      // ADD to V diffusion contribution
      const auto fd = m_stencilVals(2) - two*m_stencilVals(1) + m_stencilVals(0);
      V(smPt) += dxInvSq*m_diffusionCoeff*fd;

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (m_probEn == ::pressiodemoapps::DiffusionReaction1d::ProblemA &&
	  !m_onlyComputeVelocity)
      {
	// note that the 3 here is because of the contribution from ghost cells
	auto jvalue = -three*diffDxInvSq + twoReacCoeff*m_stencilVals(1);
	m_tripletList.push_back( Tr(smPt, uIndex, jvalue) );

	if (uIndexLeft != -1){
	  m_tripletList.push_back( Tr(smPt, uIndexLeft, diffDxInvSq) );
	}

	if (uIndexRight != -1){
	  m_tripletList.push_back( Tr(smPt, uIndexRight, diffDxInvSq) );
	}
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
    const auto & graph  = m_meshObj.graph();
    const auto & x      = m_meshObj.viewX();
    constexpr auto two  = static_cast<scalar_t>(2);

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
    const auto twoReacCoeff = m_reactionCoeff*two;
    const auto diffDxInvSq  = m_diffusionCoeff*dxInvSq;
#endif

    const auto & rows   = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (std::size_t it=0; it<rows.size(); ++it)
    {
      const auto smPt        = rows[it];
      const auto uIndex	     = graph(smPt, 0);
      const auto uIndexLeft  = graph(smPt, 1);
      const auto uIndexRight = graph(smPt, 2);

      // compute source, store into V
      m_sourceFunctor(x(uIndex), t, V(smPt));

      // ADD to V reaction contribution
      V(smPt) += m_reactionCoeff*U(uIndex)*U(uIndex);

      // ADD to V diffusion contribution
      const auto fd = U(uIndexRight) - two*U(uIndex) + U(uIndexLeft);
      V(smPt) += dxInvSq*m_diffusionCoeff*fd;

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
      if (m_probEn == ::pressiodemoapps::DiffusionReaction1d::ProblemA &&
	  !m_onlyComputeVelocity)
      {
	const auto jvalueself = -two*diffDxInvSq + twoReacCoeff*U(uIndex);
	m_tripletList.push_back( Tr(smPt, uIndex,      jvalueself) );
	m_tripletList.push_back( Tr(smPt, uIndexLeft,  diffDxInvSq) );
	m_tripletList.push_back( Tr(smPt, uIndexRight, diffDxInvSq) );
      }
#endif
    }
  }

private:
  void allocateGhosts()
  {
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    const auto numGhostValues = numDofPerCell*((stencilSize-1)/2);
    ::pressiodemoapps::resize(m_ghostLeft,  numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight, numGhostValues);
  }

private:
  const mesh_t & m_meshObj;
  int m_stencilSize = {};
  ::pressiodemoapps::DiffusionReaction1d m_probEn;
  ::pressiodemoapps::ViscousFluxReconstruction m_recEn;

  std::function<void(const scalar_t &, const scalar_t &, scalar_t &)> m_sourceFunctor;
  scalar_t m_diffusionCoeff = {};
  scalar_t m_reactionCoeff = {};

  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};
  mutable stencil_values_t m_stencilVals;
  mutable ghost_t m_ghostLeft;
  mutable ghost_t m_ghostRight;

#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
  mutable jacobian_type m_jacobian = {};
  mutable std::vector<Tr> m_tripletList;
  bool m_onlyComputeVelocity = false;
#endif
};

}}
#endif
