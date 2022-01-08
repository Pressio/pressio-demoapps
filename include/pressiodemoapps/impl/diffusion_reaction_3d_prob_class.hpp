
#ifndef PRESSIODEMOAPPS_DIFF_REACTION_3D_IMPL_HPP_
#define PRESSIODEMOAPPS_DIFF_REACTION_3D_IMPL_HPP_

#include "functor_fill_stencil.hpp"

#ifdef PRESSIODEMOAPPS_ENABLE_OPENMP
#include <omp.h>
#endif

namespace pressiodemoapps{ namespace impldiffreac{

template<class MeshType>
class EigenDiffReac3dApp
{

public:
  using index_t	      = typename MeshType::index_t;
  using scalar_type   = typename MeshType::scalar_t;
  using state_type    = Eigen::Matrix<scalar_type,Eigen::Dynamic,1>;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, index_t>;

  static constexpr int dimensionality{3};

public:
  EigenDiffReac3dApp() = delete;

  EigenDiffReac3dApp(const MeshType & meshObj,
		     ::pressiodemoapps::DiffusionReaction3d probEnum,
		     ::pressiodemoapps::ViscousFluxReconstruction recEnum,
		     scalar_type diffusion_u,
		     scalar_type diffusion_v,
		     scalar_type feedRate,
		     scalar_type killRate)
    : m_meshObj(meshObj), m_probEn(probEnum), m_recEn(recEnum)
  {
    if (m_probEn == ::pressiodemoapps::DiffusionReaction3d::GrayScott){
      setupForGrayScott(diffusion_u, diffusion_v, feedRate, killRate);
    }
  }

  EigenDiffReac3dApp(const MeshType & meshObj,
		     ::pressiodemoapps::DiffusionReaction3d probEnum,
		     ::pressiodemoapps::ViscousFluxReconstruction recEnum)
    : m_meshObj(meshObj), m_probEn(probEnum), m_recEn(recEnum)
  {
    if (m_probEn == ::pressiodemoapps::DiffusionReaction3d::GrayScott){
      setupForGrayScott(0.0002, 0.0002/2., 0.035, 0.062);
    }
  }

  state_type initialCondition() const
  {
    state_type ic(m_numDofStencilMesh);

    if (m_probEn == ::pressiodemoapps::DiffusionReaction3d::GrayScott){
      return initialConditionGS1();
    }

    return ic;
  }

protected:
  state_type initialConditionGS1() const
  {
    state_type ic(m_numDofStencilMesh);

    const auto &x= m_meshObj.viewX();
    const auto &y= m_meshObj.viewY();
    const auto &z= m_meshObj.viewZ();

    for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
      {
	const auto ind = i*m_numDofPerCell;
	if (std::abs(x(i)) < 0.3 &&
	    std::abs(y(i)) < 0.3 &&
	    std::abs(z(i)) < 0.3)
	  {
	    ic(ind)   = 1.;  // u
	    ic(ind+1) = 1.;  // v
	  }
	else{
	  ic(ind)   = 1.; // u
	  ic(ind+1) = 0.; // v
	}

      }

    return ic;
  }

  state_type initialConditionGS2() const
  {
    state_type ic(m_numDofStencilMesh);

    const auto &x= m_meshObj.viewX();
    const auto &y= m_meshObj.viewY();
    const auto &z= m_meshObj.viewZ();

    for (int i=0; i<::pressiodemoapps::extent(x,0); ++i)
      {
	const auto ind = i*m_numDofPerCell;
	const auto distX = x(i)-0.5;
	const auto distY = y(i)-0.5;
	const auto distZ = z(i)-0.5;
	const auto xsq = distX*distX;
	const auto ysq = distY*distY;
	const auto zsq = distZ*distZ;
	const auto myR = std::sqrt(xsq+ysq+zsq);
	if (myR < 0.2)
	  {
	    ic(ind)   = 1.; // v
	    ic(ind+1) = 1.;  // v
	  }
	else{
	  ic(ind)   = 1.; // u
	  ic(ind+1) = 0.; // v
	}
      }

    return ic;
  }

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

	for (index_t i=1; i<=6; ++i){
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

    int nonZerosCountBeforeComputing = 0;
    if (J){
      nonZerosCountBeforeComputing = J->nonZeros();
    }

    if (m_probEn == ::pressiodemoapps::DiffusionReaction3d::GrayScott){
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

  template<class U_t, class V_t>
  void velocityAndOptionalJacobianGrayScott(const U_t & U,
					    const scalar_type currentTime,
					    V_t & V,
					    jacobian_type * J) const
  {
    const auto & x = m_meshObj.viewX();
    const auto & y = m_meshObj.viewY();
    const auto & z = m_meshObj.viewZ();

    constexpr auto one  = static_cast<scalar_type>(1);
    constexpr auto two  = static_cast<scalar_type>(2);

    const auto dxInvSq  = m_meshObj.dxInv()*m_meshObj.dxInv();
    const auto dyInvSq  = m_meshObj.dyInv()*m_meshObj.dyInv();
    const auto dzInvSq  = m_meshObj.dzInv()*m_meshObj.dzInv();
    const auto u_diffDxInvSq  = m_diffusionCoeff_u*dxInvSq;
    const auto u_diffDyInvSq  = m_diffusionCoeff_u*dyInvSq;
    const auto u_diffDzInvSq  = m_diffusionCoeff_u*dzInvSq;

    const auto v_diffDxInvSq  = m_diffusionCoeff_v*dxInvSq;
    const auto v_diffDyInvSq  = m_diffusionCoeff_v*dyInvSq;
    const auto v_diffDzInvSq  = m_diffusionCoeff_v*dzInvSq;

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
      const auto stateIndexBottom = graph(smPt, 5)*m_numDofPerCell;
      const auto stateIndexTop    = graph(smPt, 6)*m_numDofPerCell;

      const auto thisCell_u = U(stateIndex);
      const auto thisCell_v = U(stateIndex+1);
      const auto uvSquared = thisCell_u * thisCell_v * thisCell_v;

      V(vIndexCurrentCellFirstDof) = m_feedRate * (one - thisCell_u)
	- uvSquared
	+ u_diffDxInvSq*( U(stateIndexRight) - two*U(stateIndex) + U(stateIndexLeft) )
	+ u_diffDyInvSq*( U(stateIndexBack)  - two*U(stateIndex) + U(stateIndexFront) )
	+ u_diffDzInvSq*( U(stateIndexBottom)- two*U(stateIndex) + U(stateIndexTop) );

      V(vIndexCurrentCellFirstDof+1) = -(m_feedRate+ m_killRate) * thisCell_v
	+ uvSquared
	+ v_diffDxInvSq*( U(stateIndexRight+1) - two*U(stateIndex+1) + U(stateIndexLeft+1) )
	+ v_diffDyInvSq*( U(stateIndexBack+1)  - two*U(stateIndex+1) + U(stateIndexFront+1) )
	+ v_diffDzInvSq*( U(stateIndexBottom+1)- two*U(stateIndex+1) + U(stateIndexTop+1) );

      if(J){
	// \partial f_1/\partial u
	(*J).coeffRef(vIndexCurrentCellFirstDof, stateIndex) +=
	  -two*u_diffDxInvSq - two*u_diffDyInvSq - thisCell_v * thisCell_v - m_feedRate;

	// \partial f_1/\partial v
	(*J).coeffRef(vIndexCurrentCellFirstDof, stateIndex+1) -= two*thisCell_u * thisCell_v;

	(*J).coeffRef(vIndexCurrentCellFirstDof, stateIndexLeft)   += u_diffDxInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof, stateIndexRight)  += u_diffDxInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof, stateIndexFront)  += u_diffDyInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof, stateIndexBack)   += u_diffDyInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof, stateIndexBottom) += u_diffDzInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof, stateIndexTop)    += u_diffDzInvSq;

	// \partial f_2/\partial u
	(*J).coeffRef(vIndexCurrentCellFirstDof+1, stateIndex) += thisCell_v * thisCell_v;

	// \partial f_2/\partial v
	(*J).coeffRef(vIndexCurrentCellFirstDof+1, stateIndex+1) +=
	  -two*v_diffDxInvSq -two*v_diffDyInvSq + two * thisCell_u * thisCell_v - (m_feedRate+m_killRate);

	(*J).coeffRef(vIndexCurrentCellFirstDof+1, stateIndexLeft+1)  += v_diffDxInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof+1, stateIndexRight+1) += v_diffDxInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof+1, stateIndexFront+1) += v_diffDyInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof+1, stateIndexBack+1)  += v_diffDyInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof+1, stateIndexBottom+1) += v_diffDzInvSq;
	(*J).coeffRef(vIndexCurrentCellFirstDof+1, stateIndexTop+1)    += v_diffDzInvSq;
      }
    }
  }

protected:
  const MeshType & m_meshObj;
  ::pressiodemoapps::DiffusionReaction3d m_probEn;
  ::pressiodemoapps::ViscousFluxReconstruction m_recEn;

  int m_numDofPerCell = {};
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  scalar_type m_diffusionCoeff_u = {};
  scalar_type m_diffusionCoeff_v = {};
  scalar_type m_feedRate = {};
  scalar_type m_killRate = {};
};

template<class MeshType> constexpr int EigenDiffReac3dApp<MeshType>::dimensionality;

}}
#endif
