
#ifndef PRESSIODEMOAPPS_DIFF_REACTION_2D_IMPL_HPP_
#define PRESSIODEMOAPPS_DIFF_REACTION_2D_IMPL_HPP_

#include "functor_fill_stencil.hpp"
#include "diffusion_reaction_2d_ghost_filler.hpp"

namespace pressiodemoapps{ namespace impldiffreac{

// this is the default source functor
template<class scalar_type>
struct DefaultSourceF2d
{
  void operator()(const scalar_type & x,
		  const scalar_type & y,
		  const scalar_type & evaltime,
		  scalar_type & value)
  {
    (void) evaltime;
    value = std::sin(M_PI*x*(y-0.2)) * 4.*std::sin(4.*M_PI*y*x);
  }
};

template<class MeshType>
class EigenDiffReac2dApp
{

public:
  using index_t	      = typename MeshType::index_t;
  using scalar_type   = typename MeshType::scalar_t;
  using state_type    = Eigen::Matrix<scalar_type,Eigen::Dynamic,1>;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, index_t>;

  static constexpr int	   dimensionality{2};
  static constexpr index_t numDofPerCell{1};

private:
  using ghost_container_type   = Eigen::Matrix<scalar_type,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor>;
  using stencil_container_type = Eigen::Matrix<scalar_type,Eigen::Dynamic,1>;

public:
  template<class SourceT>
  EigenDiffReac2dApp(const MeshType & meshObj,
		     ::pressiodemoapps::DiffusionReaction2d probEnum,
		     ::pressiodemoapps::ViscousFluxReconstruction recEnum,
		     SourceT sf,
		     scalar_type diffusionCoeff,
		     scalar_type reactionCoeff)
    : m_meshObj(meshObj),
      m_numDofStencilMesh(m_meshObj.stencilMeshSize()*numDofPerCell),
      m_numDofSampleMesh(m_meshObj.sampleMeshSize()*numDofPerCell),
      m_probEn(probEnum), m_recEn(recEnum)
  {
    if (m_meshObj.stencilSize() != 3){
      throw std::runtime_error("DiffusionReaction2d currently, only supports 3-pt stencil");
    }

    m_sourceFunctor = sf;
    m_diffusionCoeff = diffusionCoeff;
    m_reactionCoeff = reactionCoeff;

    const auto stencilSize = reconstructionTypeToStencilSize(recEnum);
    ::pressiodemoapps::resize(m_stencilVals, numDofPerCell*stencilSize);
    allocateGhosts();
  }

  EigenDiffReac2dApp(const MeshType & meshObj,
		     ::pressiodemoapps::DiffusionReaction2d probEnum,
		     ::pressiodemoapps::ViscousFluxReconstruction recEnum,
		     scalar_type diffusionCoeff,
		     scalar_type reactionCoeff)
    : EigenDiffReac2dApp(meshObj, probEnum, recEnum,
			 DefaultSourceF2d<scalar_type>(),
			 diffusionCoeff, reactionCoeff)
  {}

  EigenDiffReac2dApp(const MeshType & meshObj,
		     ::pressiodemoapps::DiffusionReaction2d probEnum,
		     ::pressiodemoapps::ViscousFluxReconstruction recEnum)
    : EigenDiffReac2dApp(meshObj, probEnum, recEnum,
			 DefaultSourceF2d<scalar_type>(),
			 0.01, 0.01)
  {}

  state_type initialCondition() const
  {
    state_type ic(m_numDofStencilMesh);
    if (m_probEn == ::pressiodemoapps::DiffusionReaction2d::ProblemA){
      for (int i=0; i<::pressiodemoapps::extent(ic,0); ++i){
	ic(i) = scalar_type(0);
      }
    }
    return ic;
  }

protected:
  void initializeJacobian(jacobian_type & J)
  {
    J.resize(m_numDofSampleMesh, m_numDofStencilMesh);

    using Tr = Eigen::Triplet<scalar_type>;
    std::vector<Tr> trList;

    const scalar_type val0 = 0;
    const auto & graph = m_meshObj.graph();
    for (int cell=0; cell<m_meshObj.sampleMeshSize(); ++cell)
      {
	const auto ci   = graph(cell, 0)*numDofPerCell;
	trList.push_back( Tr(cell, ci, val0) );

	for (index_t i=1; i<=4; ++i){
	  const auto neighID = graph(cell, i);
	  if( neighID != -1){
	    trList.push_back( Tr(cell, neighID*numDofPerCell, val0) );
	  }
	}
      }

    J.setFromTriplets(trList.begin(), trList.end());

    // compress to make it a real Crs matrix
    if (!J.isCompressed()){
      J.makeCompressed();
    }
  }

  // note that here we MUST use a template because when doing
  // bindings, this gets deduced to be a Ref
  template<class U_t, class V_t>
  void velocityAndOptionalJacobian(const U_t & state,
				   const scalar_type currentTime,
				   V_t & V,
				   jacobian_type * J) const
  {
    fillGhostsIfNeeded(state);

    int nonZerosCountBeforeComputing = 0;
    if (J){
      nonZerosCountBeforeComputing = J->nonZeros();
      ::pressiodemoapps::set_zero(*J);
    }

    velocityAndOptionalJacobianNearBd(state, currentTime, V, J);
    velocityAndOptionalJacobianInnerCells(state, currentTime, V, J);

    if (J){
      assert(nonZerosCountBeforeComputing == J->nonZeros());
    }
  }

private:
  template<class U_t>
  void fillGhostsIfNeeded(const U_t & U) const
  {
    if (m_probEn == ::pressiodemoapps::DiffusionReaction2d::ProblemA)
    {
      using ghost_filler_t = GhostFillerProblemA2d<U_t, MeshType, ghost_container_type>;
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

  // note that here we MUST use a template because when doing
  // bindings, this gets deduced to be a Ref
  template<class U_t, class V_t>
  void velocityAndOptionalJacobianNearBd(const U_t & U,
					 const scalar_type currentTime,
					 V_t & V,
					 jacobian_type * J) const
  {
    // note that numDofPerCell == 1, so we omit it below
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;

    // stencil filler needed because we are doing cells near boundaries
    using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell,
      stencil_container_type, U_t, MeshType, ghost_container_type>;

    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    sfiller_t StencilFillerX(stencilSize, U, m_meshObj,
			     m_ghostLeft, m_ghostRight,
			     m_stencilVals, xAxis);
    sfiller_t StencilFillerY(stencilSize, U, m_meshObj,
			     m_ghostBack, m_ghostFront,
			     m_stencilVals, yAxis);

    const auto & x	    = m_meshObj.viewX();
    const auto & y          = m_meshObj.viewY();
    const auto & graph      = m_meshObj.graph();
    constexpr auto two      = static_cast<scalar_type>(2);
    constexpr auto three    = static_cast<scalar_type>(3);
    const auto dxInvSq	    = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto dyInvSq	    = m_meshObj.dyInv()*m_meshObj.dyInv();
    const auto twoReacCoeff = m_reactionCoeff*two;
    const auto diffDxInvSq  = m_diffusionCoeff*dxInvSq;
    const auto diffDyInvSq  = m_diffusionCoeff*dyInvSq;

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
      m_sourceFunctor(x(uIndex), y(uIndex), currentTime, V(smPt));

      // add to V reaction contribution
      V(smPt) += m_reactionCoeff*U(uIndex)*U(uIndex);

      // *** add X contribution of diffusion ***
      StencilFillerX(smPt, it);
      V(smPt) += dxInvSq*m_diffusionCoeff*( m_stencilVals(2)
					    -two*m_stencilVals(1)
					    +m_stencilVals(0) );

      // *** add Y contribution of diffusion ***
      StencilFillerY(smPt, it);
      V(smPt) += dyInvSq*m_diffusionCoeff*( m_stencilVals(2)
					    -two*m_stencilVals(1)
					    +m_stencilVals(0) );

      if (m_probEn == ::pressiodemoapps::DiffusionReaction2d::ProblemA && J)
      {
	auto selfValue = -two*diffDxInvSq -two*diffDyInvSq + twoReacCoeff*m_stencilVals(1);

	if (uIndexLeft != -1){
	  J->coeffRef(smPt, uIndexLeft) += diffDxInvSq;
	}else{
	  selfValue += -diffDxInvSq;
	}

	if (uIndexFront != -1){
	  J->coeffRef(smPt, uIndexFront) += diffDyInvSq;
	}else{
	  selfValue += -diffDyInvSq;
	}

	if (uIndexRight != -1){
	  J->coeffRef(smPt, uIndexRight) += diffDxInvSq;
	}
	else{
	  selfValue += -diffDxInvSq;
	}

	if (uIndexBack != -1){
	  J->coeffRef(smPt, uIndexBack) += diffDyInvSq;
	}
	else{
	  selfValue += -diffDyInvSq;
	}
	J->coeffRef(smPt, uIndex) = selfValue;
      }
    }
  }

  template<class U_t, class V_t>
  void velocityAndOptionalJacobianInnerCells(const U_t & U,
					     const scalar_type currentTime,
					     V_t & V,
					     jacobian_type * J) const
  {
    // note that numDofPerCell == 1, so we omit it below

    const auto & x      = m_meshObj.viewX();
    const auto & y      = m_meshObj.viewY();
    const auto & graph  = m_meshObj.graph();
    constexpr auto two  = static_cast<scalar_type>(2);
    const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto dyInvSq  = m_meshObj.dyInv()*m_meshObj.dyInv();
    const auto twoReacCoeff = m_reactionCoeff*two;
    const auto diffDxInvSq  = m_diffusionCoeff*dxInvSq;
    const auto diffDyInvSq  = m_diffusionCoeff*dyInvSq;

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
      m_sourceFunctor(x(uIndex), y(uIndex), currentTime, V(smPt));

      // ADD to V reaction contribution
      V(smPt) += m_reactionCoeff*U(uIndex)*U(uIndex);

      // ADD to V x diffusion contribution
      V(smPt) += dxInvSq*m_diffusionCoeff*( U(uIndexRight) - two*U(uIndex) + U(uIndexLeft) );

      // ADD to V x diffusion contribution
      V(smPt) += dyInvSq*m_diffusionCoeff*( U(uIndexFront) - two*U(uIndex) + U(uIndexBack) );

      if (m_probEn == ::pressiodemoapps::DiffusionReaction2d::ProblemA && J)
      {
	const auto jvalueself = -two*diffDxInvSq -two*diffDyInvSq + twoReacCoeff*U(uIndex);
	(*J).coeffRef(smPt, uIndex)      += jvalueself;
	(*J).coeffRef(smPt, uIndexLeft)  += diffDxInvSq;
	(*J).coeffRef(smPt, uIndexFront) += diffDyInvSq;
	(*J).coeffRef(smPt, uIndexRight) += diffDxInvSq;
	(*J).coeffRef(smPt, uIndexBack)  += diffDyInvSq;
      }
    }
  }

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

protected:
  const std::array<scalar_type, 2> normalX_{1, 0};
  const std::array<scalar_type, 2> normalY_{0, 1};

  const MeshType & m_meshObj;
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  ::pressiodemoapps::DiffusionReaction2d m_probEn;
  ::pressiodemoapps::ViscousFluxReconstruction m_recEn;

  std::function<void(const scalar_type & /*x*/,
		     const scalar_type & /*y*/,
		     const scalar_type & /*time*/,
		     scalar_type &)> m_sourceFunctor;
  scalar_type m_diffusionCoeff = {};
  scalar_type m_reactionCoeff = {};

  mutable stencil_container_type m_stencilVals;
  mutable ghost_container_type m_ghostLeft;
  mutable ghost_container_type m_ghostFront;
  mutable ghost_container_type m_ghostRight;
  mutable ghost_container_type m_ghostBack;
};

}}
#endif
