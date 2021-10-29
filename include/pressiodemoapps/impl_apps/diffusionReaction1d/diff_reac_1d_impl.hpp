
#ifndef PRESSIODEMOAPPS_DIFF_REACTION_1D_IMPL_HPP_
#define PRESSIODEMOAPPS_DIFF_REACTION_1D_IMPL_HPP_

#include "../stencil_filler.hpp"
#include "./ghost_filler.hpp"

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
  class jacobian_t
  >
class DiffReac1dAppT
{
  static_assert(std::is_same<scalar_t, typename mesh_t::scalar_t>::value, "");

public:
  using index_t		 = typename mesh_t::index_t;
  using stencil_values_t = state_t;
  using scalar_type	 = scalar_t;
  using state_type	 = state_t;
  using velocity_type	 = state_t;
  using jacobian_type	 = jacobian_t;
  using ghost_type	 = ghost_t;
#ifdef PRESSIODEMOAPPS_ENABLE_TPL_EIGEN
  typedef Eigen::Triplet<scalar_type> Tr;
#endif

  static constexpr int dimensionality{1};
  static constexpr int numDofPerCell{1};

public:

  template<class SourceT>
  DiffReac1dAppT(const mesh_t & meshObj,
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

  DiffReac1dAppT(const mesh_t & meshObj,
		 ::pressiodemoapps::DiffusionReaction1d probEnum,
		 ::pressiodemoapps::ViscousFluxReconstruction recEnum,
		 scalar_t diffusionCoeff,
		 scalar_t reactionCoeff)
    : DiffReac1dAppT(meshObj, probEnum, recEnum,
		     DefaultSourceF<scalar_t>(), diffusionCoeff, reactionCoeff)
  {}

  DiffReac1dAppT(const mesh_t & meshObj,
		 ::pressiodemoapps::DiffusionReaction1d probEnum,
		 ::pressiodemoapps::ViscousFluxReconstruction recEnum)
    : DiffReac1dAppT(meshObj, probEnum, recEnum,
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

  jacobian_type createJacobian() const{
    jacobian_type JJ(m_numDofSampleMesh, m_numDofStencilMesh);
    return JJ;
  }

  void velocity(const state_type & state,
		const scalar_type timeValue,
		velocity_type & velo) const
  {
    velocityImpl(state, timeValue, velo);
  }

  void jacobian(const state_type & state,
		const scalar_type time,
		jacobian_type & jacobian) const
  {
    jacobianImpl(state, time, jacobian);
  }

private:
  void velocityImpl(const state_type & U,
		    const scalar_type t,
		    velocity_type & V) const
  {
    // note that numDofPerCell == 1, so we omit it below

    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);

    // only need ghosts for specific problems
    if (m_probEn == ::pressiodemoapps::DiffusionReaction1d::ProblemA)
    {
      using ghost_filler_t  = ::pressiodemoapps::impldiffreac::GhostFillerHomDirich1d<
	numDofPerCell, state_type, mesh_t, ghost_t>;
      ghost_filler_t ghF(stencilSize, U, m_meshObj, m_ghostLeft, m_ghostRight);
      ghF();
    }

    // stencil filler always needed
    using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_values_t, state_type, mesh_t, ghost_t>;
    sfiller_t StencilFiller(stencilSize, U, m_meshObj,
			    m_ghostLeft, m_ghostRight,
			    m_stencilVals);

    const auto dxInvSq = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    const auto & graph = m_meshObj.graph();
    const auto & x = m_meshObj.viewX();

    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
    {
      StencilFiller(smPt);
      const auto uIndex = graph(smPt, 0);
      const auto rhsIndex = smPt;

      // compute source
      m_sourceFunctor(x(uIndex), t, V(rhsIndex));

      // add reaction
      V(rhsIndex) += m_reactionCoeff*U(uIndex)*U(uIndex);

      constexpr auto two = static_cast<scalar_t>(2);
      const auto fd = m_stencilVals(2) - two*m_stencilVals(1) + m_stencilVals(0);
      V(rhsIndex) += dxInvSq*m_diffusionCoeff*fd;
    }
  }

  void jacobianImpl(const state_type & state,
		    const scalar_type time,
		    jacobian_type & jacobian) const
  {
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);

    using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_values_t, state_type, mesh_t, ghost_t>;
    sfiller_t StencilFiller(stencilSize, state, m_meshObj,
			    m_ghostLeft, m_ghostRight,
			    m_stencilVals);

    const auto dxInvSq = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    const auto & graph = m_meshObj.graph();
    const auto & x = m_meshObj.viewX();

    constexpr auto two   = static_cast<scalar_t>(2);
    constexpr auto three = static_cast<scalar_t>(3);
    const auto twoReacCoeff = m_reactionCoeff*two;
    const auto diffDxInvSq  = m_diffusionCoeff*dxInvSq;

    m_tripletList.clear();
    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
    {
      StencilFiller(smPt);
      const auto rowIndex    = smPt;
      const auto uIndex      = graph(smPt, 0);
      const auto uIndexLeft  = graph(smPt, 1);
      const auto uIndexRight = graph(smPt, 2);

      auto jvalue = -two*diffDxInvSq + twoReacCoeff*m_stencilVals(1);

      if (uIndexLeft != -1){
	jvalue += -diffDxInvSq;
	m_tripletList.push_back( Tr(rowIndex, uIndexLeft, diffDxInvSq) );
      }

      if (uIndexRight != -1){
	jvalue += -diffDxInvSq;
	m_tripletList.push_back( Tr(rowIndex, uIndexRight, diffDxInvSq) );
      }

      m_tripletList.push_back( Tr(rowIndex, uIndex, jvalue) );
    }

    jacobian.setFromTriplets(m_tripletList.begin(), m_tripletList.end());
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
  mutable std::vector<Tr> m_tripletList;
#endif
};

template<class T1, class T2, class T3, class T4, class T5, class T6>
constexpr int DiffReac1dAppT<T1, T2, T3, T4, T5, T6>::dimensionality;

template<class T1, class T2, class T3, class T4, class T5, class T6>
constexpr int DiffReac1dAppT<T1, T2, T3, T4, T5, T6>::numDofPerCell;

}}
#endif
