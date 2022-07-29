
#ifndef PRESSIODEMOAPPS_DIFF_REACTION_1D_IMPL_HPP_
#define PRESSIODEMOAPPS_DIFF_REACTION_1D_IMPL_HPP_

#include "functor_fill_stencil.hpp"
#include "diffusion_reaction_1d_ghost_filler.hpp"

namespace pressiodemoapps{
namespace impldiffusionreaction1d{

// tags are used inside he public create function: create_problem_...()
// to dispatch to the proper problem
// so add new ones if a new problem is added
struct TagProblemA{};

template<class MeshType>
class EigenApp
{

public:
  // required public aliases
  using index_t	      = typename MeshType::index_t;
  using scalar_type   = typename MeshType::scalar_t;
  using state_type    = Eigen::Matrix<scalar_type,Eigen::Dynamic,1>;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, index_t>;

private:
  static constexpr int dimensionality{1};
  using stencil_values_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using ghost_container_type = Eigen::Matrix<scalar_type,
					     Eigen::Dynamic,
					     Eigen::Dynamic,
					     Eigen::RowMajor>;

public:
  template<class SourceT>
  EigenApp(TagProblemA /*tag*/,
	   const MeshType & meshObj,
	   ::pressiodemoapps::ViscousFluxReconstruction recEnum,
	   SourceT sourceFunctor,
	   scalar_type diffusionCoeff,
	   scalar_type reactionCoeff)
    : m_numDofPerCell(1),
      m_probEn(::pressiodemoapps::DiffusionReaction1d::ProblemA),
      m_viscousFluxRecEn(recEnum),
      m_meshObj(meshObj),
      m_numDofStencilMesh(m_meshObj.stencilMeshSize()),
      m_numDofSampleMesh(m_meshObj.sampleMeshSize()),
      m_probA_diffusionCoeff(diffusionCoeff),
      m_probA_reactionCoeff(reactionCoeff)
  {

    if (m_meshObj.stencilSize() != 3){
      throw std::runtime_error("DiffusionReaction1d currently only supports 3-pt stencil");
    }

    m_probA_sourceFunctor = sourceFunctor;
    allocateGhosts();
  }

  state_type initialCondition() const
  {
    state_type ic(m_numDofStencilMesh);
    if (m_probEn == ::pressiodemoapps::DiffusionReaction1d::ProblemA){
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

    constexpr auto val0 = static_cast<scalar_type>(0);
    const auto & graph = m_meshObj.graph();
    for (int cell=0; cell<m_meshObj.sampleMeshSize(); ++cell)
      {
	const auto ci   = graph(cell, 0)*m_numDofPerCell;
	trList.push_back( Tr(cell, ci, val0) );

	// this is 1d and for stencil = 3, we only need left and right neighbors
	for (index_t i=1; i<=2; ++i){
	  const auto neighID = graph(cell, i);
	  if( neighID != -1){
	    trList.push_back( Tr(cell, neighID*m_numDofPerCell, val0) );
	  }
	}
      }

    J.setFromTriplets(trList.begin(), trList.end());

    // compress to make it a real Crs matrix
    if (!J.isCompressed()){
      J.makeCompressed();
    }
  }

  template<class U_t, class V_t>
  void velocityAndOptionalJacobian(const U_t & state,
				   const scalar_type currentTime,
				   V_t & V,
				   jacobian_type * J) const
  {
    fillGhostsIfNeeded(state, currentTime);

    int nonZerosCountBeforeComputing = 0;
    if (J){
      nonZerosCountBeforeComputing = J->nonZeros();
      ::pressiodemoapps::set_zero(*J);
    }

    velocityAndOptionalJacobianNearBd(state, currentTime, V, J);
    velocityAndOptionalJacobianInnerCells(state, currentTime, V, J);

    if (J){
      assert(nonZerosCountBeforeComputing == J->nonZeros());
      (void) nonZerosCountBeforeComputing;
    }
  }

private:
  template<class U_t>
  void fillGhostsIfNeeded(const U_t & U,
			  scalar_type /*currTime*/) const
  {
    if (m_probEn == ::pressiodemoapps::DiffusionReaction1d::ProblemA){
      using ghost_filler_t  = GhostFillerProblemA1d<U_t, MeshType, ghost_container_type>;
      const auto stencilSizeNeeded = reconstructionTypeToStencilSize(m_viscousFluxRecEn);
      ghost_filler_t ghF(stencilSizeNeeded, U, m_meshObj, m_ghostLeft, m_ghostRight);

      const auto & rowsBd = m_meshObj.graphRowsOfCellsNearBd();
      for (decltype(rowsBd.size()) it=0; it<rowsBd.size(); ++it){
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
    // but will need to change if we implement another problem


    const auto stencilSize = reconstructionTypeToStencilSize(m_viscousFluxRecEn);
    stencil_values_t stencilVals(m_numDofPerCell*stencilSize);

    // stencil filler needed because we are doing cells near boundaries
    using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, stencil_values_t, U_t, MeshType, ghost_container_type>;
    sfiller_t StencilFiller(stencilSize, U, m_meshObj,
			    m_ghostLeft, m_ghostRight, stencilVals);

    const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto & graph  = m_meshObj.graph();
    const auto & x      = m_meshObj.viewX();
    constexpr auto two  = static_cast<scalar_type>(2);
    constexpr auto three= static_cast<scalar_type>(3);
    const auto twoReacCoeff = m_probA_reactionCoeff*two;
    const auto diffDxInvSq  = m_probA_diffusionCoeff*dxInvSq;

    const auto & rows   = m_meshObj.graphRowsOfCellsNearBd();
    for (std::size_t it=0; it<rows.size(); ++it)
    {
      const auto smPt        = rows[it];
      const auto uIndex	     = graph(smPt, 0);
      const auto uIndexLeft  = graph(smPt, 1);
      const auto uIndexRight = graph(smPt, 2);

      StencilFiller(smPt, it, 1 /*1 dof/cell */);

      // compute source, put contribution into V
      m_probA_sourceFunctor(x(uIndex), currentTime, V(smPt));

      // ADD to V reaction contribution
      V(smPt) += m_probA_reactionCoeff*U(uIndex)*U(uIndex);

      // ADD to V diffusion contribution
      const auto fd = stencilVals(2) - two*stencilVals(1) + stencilVals(0);
      V(smPt) += dxInvSq*m_probA_diffusionCoeff*fd;

      if (m_probEn == ::pressiodemoapps::DiffusionReaction1d::ProblemA && J)
      {
	// the 3 here stems from the contribution from ghost cell
	auto jvalue = -three*diffDxInvSq + twoReacCoeff*stencilVals(1);
	(*J).coeffRef(smPt, uIndex) += jvalue;

	if (uIndexLeft != -1){
	  (*J).coeffRef(smPt, uIndexLeft) += diffDxInvSq;
	}

	if (uIndexRight != -1){
	  (*J).coeffRef(smPt, uIndexRight) += diffDxInvSq;
	}
      }
    }
  }

  template<class U_t, class V_t>
  void velocityAndOptionalJacobianInnerCells(const U_t & U,
					     const scalar_type currentTime,
					     V_t & V,
					     jacobian_type * J) const
  {

    const auto & x      = m_meshObj.viewX();
    const auto & graph  = m_meshObj.graph();
    constexpr auto two  = static_cast<scalar_type>(2);
    const auto twoReacCoeff = m_probA_reactionCoeff*two;
    const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto diffDxInvSq  = m_probA_diffusionCoeff*dxInvSq;

    const auto & rows   = m_meshObj.graphRowsOfCellsAwayFromBd();
    for (std::size_t it=0; it<rows.size(); ++it)
    {
      const auto smPt        = rows[it];
      const auto uIndex	     = graph(smPt, 0);
      const auto uIndexLeft  = graph(smPt, 1);
      const auto uIndexRight = graph(smPt, 2);

      // compute source, store into V
      m_probA_sourceFunctor(x(uIndex), currentTime, V(smPt));

      // ADD to V reaction contribution
      V(smPt) += m_probA_reactionCoeff*U(uIndex)*U(uIndex);

      // ADD to V diffusion contribution
      const auto fd = U(uIndexRight) - two*U(uIndex) + U(uIndexLeft);
      V(smPt) += dxInvSq*m_probA_diffusionCoeff*fd;

      if (m_probEn == ::pressiodemoapps::DiffusionReaction1d::ProblemA && J)
      {
	const auto jvalueself = -two*diffDxInvSq + twoReacCoeff*U(uIndex);
	(*J).coeffRef(smPt, uIndex)      += jvalueself;
	(*J).coeffRef(smPt, uIndexLeft)  += diffDxInvSq;
	(*J).coeffRef(smPt, uIndexRight) += diffDxInvSq;
      }
    }
  }

  void allocateGhosts()
  {
    const auto stencilSize    = reconstructionTypeToStencilSize(m_viscousFluxRecEn);
    const auto numGhostValues = m_numDofPerCell*((stencilSize-1)/2);
    const index_t s1 = m_meshObj.numCellsBd();
    ::pressiodemoapps::resize(m_ghostLeft,  s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight, s1, numGhostValues);
  }

protected:
  // common to all problems
  int m_numDofPerCell = {};
  ::pressiodemoapps::DiffusionReaction1d m_probEn;
  ::pressiodemoapps::ViscousFluxReconstruction m_viscousFluxRecEn;
  const MeshType & m_meshObj;
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};
  mutable ghost_container_type m_ghostLeft;
  mutable ghost_container_type m_ghostRight;

  // parameters specific to problems
  // will need to handle this better later
  std::function<void(const scalar_type &,
		     const scalar_type &,
		     scalar_type &)> m_probA_sourceFunctor;

  scalar_type m_probA_diffusionCoeff = {};
  scalar_type m_probA_reactionCoeff = {};
};

template<class MeshType> constexpr int EigenApp<MeshType>::dimensionality;

}}
#endif
