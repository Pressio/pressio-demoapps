
#ifndef PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_
#define PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_

// note that the code below is pretty ugly for now, but it works.
// this is inside impl namespace for a reason, and will need to be
// improved later on but we have a starting point.

namespace pressiodemoapps{ namespace ad{ namespace impl{

template<
  class scalar_t,
  class mesh_t,
  class state_t,
  class velo_t,
  class ghost_t
  >
class LinearAdvT
{
  static_assert
  (std::is_same<scalar_t, typename mesh_t::scalar_t>::value, "");

public:
  using index_t		= typename mesh_t::index_t;
  using scalar_type	= scalar_t;
  using state_type	= state_t;
  using velocity_type	= state_t;
  using ghost_type	= ghost_t;

public:
  LinearAdvT(const mesh_t & meshObj)
    : m_meshObj(meshObj)
  {
    // here we have one dof per cell
    m_numDofStencilMesh = m_meshObj.stencilMeshSize();
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize();
    m_stencilVals.resize(m_meshObj.stencilSize());
  }

  index_t totalDofSampleMesh() const{
    return m_numDofSampleMesh;
  }

  index_t totalDofStencilMesh() const{
    return m_numDofStencilMesh;
  }

  velocity_type createVelocity() const {
    velocity_type V(m_numDofSampleMesh);
    return V;
  }

  void velocity(const state_type & U,
		const scalar_type t,
		velocity_type & V) const
  {
    scalar_type FL{0};
    scalar_type FR{0};
    scalar_type uMinusHalfNeg{0};
    scalar_type uMinusHalfPos{0};
    scalar_type uPlusHalfNeg {0};
    scalar_type uPlusHalfPos {0};

    const auto dxInv = m_meshObj.dxInv();
    const auto & x = m_meshObj.viewX();
    const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    const auto & graph = m_meshObj.graph();

    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
    {
      //const auto & gigi = graph(smPt, 0);
      // const auto & gigi2 = thisCellAdList(0);
      auto thisCellAdList = pressiodemoapps::impl::neighbors(m_meshObj, smPt);
      const auto cellGID = thisCellAdList(0);

      const index_t w0 = thisCellAdList(1);
      //m_stencilVals[0] = U[w0];
      //fillStencil(smPt, cellGID, thisCellAdList, U);

    //   reconstruct(m_stencilVals,
    // 		  uMinusHalfNeg, uMinusHalfPos,
    // 		  uPlusHalfNeg, uPlusHalfPos);

    //   linAdvRusanovFlux(FL, uMinusHalfNeg, uMinusHalfPos);
    //   linAdvRusanovFlux(FR, uPlusHalfNeg,  uPlusHalfPos);

    //   V(cellGID) = dxInv*(FL - FR);
    }
  }

private:
  template<class adlist_t, class uindex_t>
  void fillStencil(index_t & smPt,
		   const uindex_t uIndex,
		   const adlist_t & adList,
		   const state_type & U) const
  {
    const index_t w0 = adList(1);
    const index_t e0 = adList(2);
    const index_t w1 = adList(3);
    const index_t e1 = adList(4);
    const index_t w2 = adList(5);
    const index_t e2 = adList(6);
    const auto stencilSize = m_meshObj.stencilSize();

    if (stencilSize ==3)
    {
      m_stencilVals[0] = U[w0];
      m_stencilVals[1] = U[uIndex];
      m_stencilVals[2] = U[e0];
    }
    else if (stencilSize ==7){
      m_stencilVals[0] = U[w2];
      m_stencilVals[1] = U[w1];
      m_stencilVals[2] = U[w0];
      m_stencilVals[3] = U[uIndex];
      m_stencilVals[4] = U[e0];
      m_stencilVals[5] = U[e1];
      m_stencilVals[6] = U[e2];
    }
  }

  template<class T, class T2>
  void reconstruct(const T2 & m_stencilVals,
		   T & uMinusHalfNeg,
		   T & uMinusHalfPos,
		   T & uPlusHalfNeg,
		   T & uPlusHalfPos) const
  {

    const auto stencilSize = m_meshObj.stencilSize();
    if (stencilSize==3)
    {
      uMinusHalfNeg = m_stencilVals[0];
      uMinusHalfPos = m_stencilVals[1];
      uPlusHalfNeg  = m_stencilVals[1];
      uPlusHalfPos  = m_stencilVals[2];
    }
    else if (stencilSize==7)
    {
      pressiodemoapps::weno5(uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg, uPlusHalfPos,
			     m_stencilVals);
    }
  }

private:
  const mesh_t & m_meshObj;
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};
  mutable std::vector<scalar_type> m_stencilVals;
};

}}}//end namespace
#endif
