
#ifndef PRESSIODEMOAPPS_DIFF_REACTION_2D_IMPL_HPP_
#define PRESSIODEMOAPPS_DIFF_REACTION_2D_IMPL_HPP_

#include "functor_fill_stencil.hpp"
#include "diffusion_reaction_2d_ghost_filler.hpp"

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#include <omp.h>
#endif

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

  static constexpr int dimensionality{2};

private:
  using ghost_container_type   = Eigen::Matrix<scalar_type,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor>;
  using stencil_container_type = Eigen::Matrix<scalar_type,Eigen::Dynamic,1>;

public:
  // constructor valid for ProblemA only
  template<class SourceT>
  EigenDiffReac2dApp(const MeshType & meshObj,
		     ::pressiodemoapps::DiffusionReaction2d probEnum,
		     ::pressiodemoapps::ViscousFluxReconstruction recEnum,
		     SourceT sf,
		     scalar_type diffusionCoeff,
		     scalar_type reactionCoeff)
    : m_meshObj(meshObj), m_probEn(probEnum), m_recEn(recEnum)
  {
    if (m_meshObj.stencilSize() != 3){
      throw std::runtime_error("DiffusionReaction2d currently, only supports 3-pt stencil");
    }

    if (m_probEn != ::pressiodemoapps::DiffusionReaction2d::ProblemA){
      throw std::runtime_error("EigenDiffReac2dApp: constructor valid for ProblemA only");
    }

    setupForProblemA(sf, diffusionCoeff, reactionCoeff);
  }

  // constructor valid for ProblemA only
  EigenDiffReac2dApp(const MeshType & meshObj,
		     ::pressiodemoapps::DiffusionReaction2d probEnum,
		     ::pressiodemoapps::ViscousFluxReconstruction recEnum,
		     scalar_type diffusionCoeff,
		     scalar_type reactionCoeff)
    : EigenDiffReac2dApp(meshObj, probEnum, recEnum,
			 DefaultSourceF2d<scalar_type>(),
			 diffusionCoeff, reactionCoeff)
  {}

  // constructor valid for GrayScott
  EigenDiffReac2dApp(const MeshType & meshObj,
		     ::pressiodemoapps::DiffusionReaction2d probEnum,
		     ::pressiodemoapps::ViscousFluxReconstruction recEnum,
		     scalar_type diffusion_u,
		     scalar_type diffusion_v,
		     scalar_type feedRate,
		     scalar_type killRate)
    : m_meshObj(meshObj), m_probEn(probEnum), m_recEn(recEnum)
  {
    if (m_probEn != ::pressiodemoapps::DiffusionReaction2d::GrayScott){
      throw std::runtime_error("EigenDiffReac2dApp: constructor valid for GrayScott only");
    }

    setupForGrayScott(diffusion_u, diffusion_v, feedRate, killRate);
  }

  // constructor valid for both ProblemA and GrayScott
  // sets all parameters to default
  EigenDiffReac2dApp(const MeshType & meshObj,
		     ::pressiodemoapps::DiffusionReaction2d probEnum,
		     ::pressiodemoapps::ViscousFluxReconstruction recEnum)
    : m_meshObj(meshObj), m_probEn(probEnum), m_recEn(recEnum)
  {
    if (m_probEn == ::pressiodemoapps::DiffusionReaction2d::ProblemA){
      setupForProblemA(DefaultSourceF2d<scalar_type>(), 0.01, 0.01);
    }
    else if (m_probEn == ::pressiodemoapps::DiffusionReaction2d::GrayScott){
      setupForGrayScott(0.0002, 0.0002/4., 0.042, 0.062);
    }
  }

  state_type initialCondition() const
  {
    state_type ic(m_numDofStencilMesh);

    if (m_probEn == ::pressiodemoapps::DiffusionReaction2d::ProblemA){
      for (int i=0; i<::pressiodemoapps::extent(ic,0); ++i){
	ic(i) = scalar_type(0);
      }
    }

    if (m_probEn == ::pressiodemoapps::DiffusionReaction2d::GrayScott)
    {
      const auto &x= m_meshObj.viewX();
      const auto &y= m_meshObj.viewY();

      for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
	{
	  const auto ind = i*m_numDofPerCell;
	  if (std::abs(x(i)) < 0.1 &&
	      std::abs(y(i)) < 0.1)
	  {
	    ic(ind)   = 0.5;   // u
	    ic(ind+1) = 0.25;  // v
	  }
	  else{
	    ic(ind)   = 1.; // u
	    ic(ind+1) = 0.; // v
	  }
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

    const scalar_type zero = {0};
    const auto & graph = m_meshObj.graph();
    for (index_t smPt=0; smPt<m_meshObj.sampleMeshSize(); ++smPt)
      {
	const auto jacRowOfCurrCellFirstDof = smPt*m_numDofPerCell;
	const auto jacColOfCurrCellFirstDof = graph(smPt, 0)*m_numDofPerCell;

	for (int k=0; k<m_numDofPerCell; ++k){
	  for (int j=0; j<m_numDofPerCell; ++j){
	    trList.push_back( Tr(jacRowOfCurrCellFirstDof+k, jacColOfCurrCellFirstDof+j, zero) );
	  }
	}

	for (index_t i=1; i<=4; ++i){
	  const auto neighID = graph(smPt, i);
	  if( neighID != -1){
	    const auto ci = neighID*m_numDofPerCell;
	    for (int k=0; k<m_numDofPerCell; ++k){
	      for (int j=0; j<m_numDofPerCell; ++j){
		trList.push_back( Tr(jacRowOfCurrCellFirstDof+k, ci+j, zero) );
	      }
	    }
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
  void velocityAndOptionalJacobian(const U_t & U,
				   const scalar_type currentTime,
				   V_t & V,
				   jacobian_type * J) const
  {

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp parallel
{
#endif

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
    ::pressiodemoapps::set_zero_omp(V);
#else
    ::pressiodemoapps::set_zero(V);
#endif

    if (J){
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
      ::pressiodemoapps::set_zero_omp(*J);
#else
      ::pressiodemoapps::set_zero(*J);
#endif
    }

    fillGhostsIfNeeded(U);

    int nonZerosCountBeforeComputing = 0;
    if (J){
      nonZerosCountBeforeComputing = J->nonZeros();
    }

    if (m_probEn == ::pressiodemoapps::DiffusionReaction2d::ProblemA){
      velocityAndOptionalJacobianNearBdProblemA(U, currentTime, V, J);
      velocityAndOptionalJacobianInnerCellsProblemA(U, currentTime, V, J);
    }
    else if (m_probEn == ::pressiodemoapps::DiffusionReaction2d::GrayScott){
      // for GrayScott we have periodic BC so no special treatment needed at boundaries
      velocityAndOptionalJacobianGrayScott(U, currentTime, V, J);
    }

    if (J){
      assert(nonZerosCountBeforeComputing == J->nonZeros());
    }

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
}//end omp parallel
#endif

  }

private:
  template<class SourceT>
  void setupForProblemA(SourceT sf,
			scalar_type diffusionCoeff,
			scalar_type reactionCoeff)
  {

    m_numDofPerCell = 1;
    m_numDofStencilMesh = m_meshObj.stencilMeshSize()*m_numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize()*m_numDofPerCell;
    m_sourceFunctor  = sf;
    m_diffusionCoeff = diffusionCoeff;
    m_reactionCoeff  = reactionCoeff;
    allocateGhosts();
  }

  void setupForGrayScott(scalar_type diffusionCoeff_u,
			 scalar_type diffusionCoeff_v,
			 scalar_type feedRate,
			 scalar_type killRate)
  {

    m_numDofPerCell = 2;
    m_numDofStencilMesh = m_meshObj.stencilMeshSize()*m_numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize()*m_numDofPerCell;
    m_diffusionCoeff_u = diffusionCoeff_u;
    m_diffusionCoeff_v = diffusionCoeff_v;
    m_feedRate = feedRate;
    m_killRate = killRate;
  }

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
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
      for (int it=0; it<rowsBd.size(); ++it){
	ghF(rowsBd[it], it);
      }
    }
  }

  // note that here we MUST use a template because when doing
  // bindings, this gets deduced to be a Ref
  template<class U_t, class V_t>
  void velocityAndOptionalJacobianNearBdProblemA(const U_t & U,
						 const scalar_type currentTime,
						 V_t & V,
						 jacobian_type * J) const
  {
    // note that m_numDofPerCell == 1, so we omit it below
    constexpr int xAxis = 1;
    constexpr int yAxis = 2;

    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    stencil_container_type stencilVals(stencilSize);

    // stencil filler needed because we are doing cells near boundaries
    using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, 1, stencil_container_type, U_t,
      MeshType, ghost_container_type>;

    sfiller_t StencilFillerX(stencilSize, U, m_meshObj,
			     m_ghostLeft, m_ghostRight,
			     stencilVals, xAxis);
    sfiller_t StencilFillerY(stencilSize, U, m_meshObj,
			     m_ghostBack, m_ghostFront,
			     stencilVals, yAxis);

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
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
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
      V(smPt) += dxInvSq*m_diffusionCoeff*( stencilVals(2)
					    -two*stencilVals(1)
					    +stencilVals(0) );

      // *** add Y contribution of diffusion ***
      StencilFillerY(smPt, it);
      V(smPt) += dyInvSq*m_diffusionCoeff*( stencilVals(2)
					    -two*stencilVals(1)
					    +stencilVals(0) );

      if (J)
      {
	auto selfValue = -two*diffDxInvSq -two*diffDyInvSq + twoReacCoeff*stencilVals(1);

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
  void velocityAndOptionalJacobianInnerCellsProblemA(const U_t & U,
						     const scalar_type currentTime,
						     V_t & V,
						     jacobian_type * J) const
  {
    // note that for ProblemA, m_numDofPerCell == 1, so we omit it below

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
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
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

      if (J){
	const auto jvalueself = -two*diffDxInvSq -two*diffDyInvSq + twoReacCoeff*U(uIndex);
	(*J).coeffRef(smPt, uIndex)      += jvalueself;
	(*J).coeffRef(smPt, uIndexLeft)  += diffDxInvSq;
	(*J).coeffRef(smPt, uIndexFront) += diffDyInvSq;
	(*J).coeffRef(smPt, uIndexRight) += diffDxInvSq;
	(*J).coeffRef(smPt, uIndexBack)  += diffDyInvSq;
      }
    }
  }

  template<class U_t, class V_t>
  void velocityAndOptionalJacobianGrayScott(const U_t & U,
					    const scalar_type currentTime,
					    V_t & V,
					    jacobian_type * J) const
  {
    const auto & x      = m_meshObj.viewX();
    const auto & y      = m_meshObj.viewY();

    constexpr auto one  = static_cast<scalar_type>(1);
    constexpr auto two  = static_cast<scalar_type>(2);

    const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto dyInvSq  = m_meshObj.dyInv()*m_meshObj.dyInv();
    const auto u_diffDxInvSq  = m_diffusionCoeff_u*dxInvSq;
    const auto u_diffDyInvSq  = m_diffusionCoeff_u*dyInvSq;
    const auto v_diffDxInvSq  = m_diffusionCoeff_v*dxInvSq;
    const auto v_diffDyInvSq  = m_diffusionCoeff_v*dyInvSq;

    const auto & graph  = m_meshObj.graph();
#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif

    for (index_t smPt=0; smPt<m_meshObj.sampleMeshSize(); ++smPt)
    {
      const auto vIndexCurrentCellFirstDof = smPt*m_numDofPerCell;

      const auto stateIndex	  = graph(smPt, 0)*m_numDofPerCell;
      const auto stateIndexLeft   = graph(smPt, 1)*m_numDofPerCell;
      const auto stateIndexFront  = graph(smPt, 2)*m_numDofPerCell;
      const auto stateIndexRight  = graph(smPt, 3)*m_numDofPerCell;
      const auto stateIndexBack   = graph(smPt, 4)*m_numDofPerCell;

      const auto thisCell_u = U(stateIndex);
      const auto thisCell_v = U(stateIndex+1);
      const auto uvSquared = thisCell_u * thisCell_v * thisCell_v;

      V(vIndexCurrentCellFirstDof) = m_feedRate * (one - thisCell_u)
	- uvSquared
	+ u_diffDxInvSq*( U(stateIndexRight) - two*U(stateIndex) + U(stateIndexLeft) )
	+ u_diffDyInvSq*( U(stateIndexBack)  - two*U(stateIndex) + U(stateIndexFront) );

      V(vIndexCurrentCellFirstDof+1) = -(m_feedRate+ m_killRate) * thisCell_v
	+ uvSquared
	+ v_diffDxInvSq*( U(stateIndexRight+1) - two*U(stateIndex+1) + U(stateIndexLeft+1) )
	+ v_diffDyInvSq*( U(stateIndexBack+1)  - two*U(stateIndex+1) + U(stateIndexFront+1) );

      if(J){
	// \partial f_1/\partial u
	(*J).coeffRef(vIndexCurrentCellFirstDof, stateIndex) +=
	  -two*u_diffDxInvSq - two*u_diffDyInvSq - thisCell_v * thisCell_v - m_feedRate;

	// \partial f_1/\partial v
	(*J).coeffRef(vIndexCurrentCellFirstDof, stateIndex+1) -= two*thisCell_u * thisCell_v;

	(*J).coeffRef(vIndexCurrentCellFirstDof, stateIndexLeft)  += u_diffDxInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof, stateIndexFront) += u_diffDyInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof, stateIndexRight) += u_diffDxInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof, stateIndexBack)  += u_diffDyInvSq;

	// \partial f_2/\partial u
	(*J).coeffRef(vIndexCurrentCellFirstDof+1, stateIndex) += thisCell_v * thisCell_v;

	// \partial f_2/\partial v
	(*J).coeffRef(vIndexCurrentCellFirstDof+1, stateIndex+1) +=
	  -two*v_diffDxInvSq -two*v_diffDyInvSq + two * thisCell_u * thisCell_v - (m_feedRate+m_killRate);

	(*J).coeffRef(vIndexCurrentCellFirstDof+1, stateIndexLeft+1)  += v_diffDxInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof+1, stateIndexFront+1) += v_diffDyInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof+1, stateIndexRight+1) += v_diffDxInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof+1, stateIndexBack+1)  += v_diffDyInvSq;
      }
    }
  }


  void allocateGhosts()
  {
    const auto stencilSize    = reconstructionTypeToStencilSize(m_recEn);
    const auto numGhostValues = m_numDofPerCell*((stencilSize-1)/2);

    const index_t s1 = m_meshObj.numCellsBd();
    ::pressiodemoapps::resize(m_ghostLeft,  s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostFront, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight, s1, numGhostValues);
    ::pressiodemoapps::resize(m_ghostBack,  s1, numGhostValues);
  }

protected:
  const MeshType & m_meshObj;
  ::pressiodemoapps::DiffusionReaction2d m_probEn;
  ::pressiodemoapps::ViscousFluxReconstruction m_recEn;

  int m_numDofPerCell         = {};
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  // members needed for problemA
  std::function<void(const scalar_type & /*x*/,
		     const scalar_type & /*y*/,
		     const scalar_type & /*time*/,
		     scalar_type &)> m_sourceFunctor;
  scalar_type m_diffusionCoeff = {};
  scalar_type m_reactionCoeff = {};
  mutable ghost_container_type m_ghostLeft;
  mutable ghost_container_type m_ghostFront;
  mutable ghost_container_type m_ghostRight;
  mutable ghost_container_type m_ghostBack;

  // members needed for Gray-Scott
  scalar_type m_diffusionCoeff_u = {};
  scalar_type m_diffusionCoeff_v = {};
  scalar_type m_feedRate = {};
  scalar_type m_killRate = {};
};

template<class MeshType> constexpr int EigenDiffReac2dApp<MeshType>::dimensionality;

}}
#endif
