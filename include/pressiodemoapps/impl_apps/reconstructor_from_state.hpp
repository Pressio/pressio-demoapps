
#ifndef PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STATE_HPP_
#define PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STATE_HPP_

namespace pressiodemoapps{ namespace impl{

template<class edge_rec_t, class state_t, class mesh_t>
struct _ReconstructorMembers
{
  _ReconstructorMembers() = delete;

  _ReconstructorMembers(::pressiodemoapps::InviscidFluxReconstruction recEn,
			const state_t & state,
			const mesh_t & meshObj,
			edge_rec_t & uMinusHalfNeg,
			edge_rec_t & uMinusHalfPos,
			edge_rec_t & uPlusHalfNeg,
			edge_rec_t & uPlusHalfPos)
    : m_recEn(recEn),
      m_state(state),
      m_meshObj(meshObj),
      m_uMinusHalfNeg(uMinusHalfNeg),
      m_uMinusHalfPos(uMinusHalfPos),
      m_uPlusHalfNeg(uPlusHalfNeg),
      m_uPlusHalfPos(uPlusHalfPos)
  {}

  _ReconstructorMembers(const int axis,
			::pressiodemoapps::InviscidFluxReconstruction recEn,
			const state_t & state,
			const mesh_t & meshObj,
			edge_rec_t & uMinusHalfNeg,
			edge_rec_t & uMinusHalfPos,
			edge_rec_t & uPlusHalfNeg,
			edge_rec_t & uPlusHalfPos)
    : m_axis(axis),
      m_recEn(recEn),
      m_state(state),
      m_meshObj(meshObj),
      m_uMinusHalfNeg(uMinusHalfNeg),
      m_uMinusHalfPos(uMinusHalfPos),
      m_uPlusHalfNeg(uPlusHalfNeg),
      m_uPlusHalfPos(uPlusHalfPos)
  {}

  int m_axis = -1;
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  edge_rec_t & m_uMinusHalfNeg;
  edge_rec_t & m_uMinusHalfPos;
  edge_rec_t & m_uPlusHalfNeg;
  edge_rec_t & m_uPlusHalfPos;

};

template<int dim, class edge_rec_t, class state_t, class mesh_t>
class ReconstructorFromState;

// specialize for 1d
template<class edge_rec_t, class state_t, class mesh_t>
class ReconstructorFromState<1, edge_rec_t, state_t, mesh_t>
{

public:
  ReconstructorFromState() = delete;

  template<typename ...Args>
  ReconstructorFromState(Args && ... args)
    : m_memb(std::forward<Args>(args)...)
  {}

  template<int ndpc>
  void operator()(int smPt)
  {
    const auto & graph  = m_memb.m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0)*ndpc;
    const auto axis = m_memb.m_axis;

    /*
      For 1d, graph is:
      ptID left0 right0 [left1 right1 left2 right2]
       0     1     2    [  3     4      5     6   ]

      where [] are optionally present, depending on the stencil chosen.
    */

    const auto l0  = graph(smPt, 1);
    const auto r0  = graph(smPt, 2);
    const auto l0i = l0*ndpc;
    const auto r0i = r0*ndpc;

    switch(m_memb.m_recEn)
      {
      case InviscidFluxReconstruction::FirstOrder:
	{
	  for (int i=0; i<ndpc; ++i){
	    m_memb.m_uMinusHalfNeg(i) = m_memb.m_state(l0i+i);
	    m_memb.m_uMinusHalfPos(i) = m_memb.m_state(uIndex+i);
	    m_memb.m_uPlusHalfNeg(i)  = m_memb.m_state(uIndex+i);
	    m_memb.m_uPlusHalfPos(i)  = m_memb.m_state(r0i+i);
	  }
	  break;
	}

      case InviscidFluxReconstruction::Weno3:
	{
	  const auto l1  = graph(smPt, 3);
	  const auto r1  = graph(smPt, 4);
	  const auto l1i = l1*ndpc;
	  const auto r1i = r1*ndpc;

	  pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg,
				 m_memb.m_uMinusHalfPos,
				 m_memb.m_uPlusHalfNeg,
				 m_memb.m_uPlusHalfPos,
				 m_memb.m_state,
				 l1i, l0i, r0i, r1i,
				 uIndex,
				 ndpc);
	  break;
	}

      case InviscidFluxReconstruction::Weno5:
	{
	  const auto l1  = graph(smPt, 3);
	  const auto r1  = graph(smPt, 4);
	  const auto l2  = graph(smPt, 5);
	  const auto r2  = graph(smPt, 6);
	  const auto l1i = l1*ndpc;
	  const auto r1i = r1*ndpc;
	  const auto l2i = l2*ndpc;
	  const auto r2i = r2*ndpc;

	  pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg,
				 m_memb.m_uMinusHalfPos,
				 m_memb.m_uPlusHalfNeg,
				 m_memb.m_uPlusHalfPos,
				 m_memb.m_state,
				 l2i, l1i, l0i,
				 r0i, r1i, r2i,
				 uIndex,
				 ndpc);
	  break;
	}
      }
  }

private:
  using members_t = _ReconstructorMembers<edge_rec_t, state_t, mesh_t>;
  members_t m_memb;
};

// specialize for 2d
template<class edge_rec_t, class state_t, class mesh_t>
class ReconstructorFromState<2, edge_rec_t, state_t, mesh_t>
{

public:
  ReconstructorFromState() = delete;

  template<typename ...Args>
  ReconstructorFromState(Args && ... args)
    : m_memb(std::forward<Args>(args)...)
  {}

  template<int ndpc>
  void operator()(int smPt)
  {
    const auto & graph  = m_memb.m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0)*ndpc;
    const auto axis = m_memb.m_axis;

    /*
      For 2d, graph is:
      ptID left0 front0 right0 back0 [left1 front1 right1 back1 left2 front2 right2 back2]
       0     1     2      3     4    [  5     6      7      8    9      10     11    12  ]

      where [] are optionally present, depending on the stencil chosen.
    */

    switch(m_memb.m_recEn)
      {
      case InviscidFluxReconstruction::FirstOrder:
	{
	  const auto l0  = (axis == 1) ? graph(smPt, 1) : graph(smPt, 4);
	  const auto r0  = (axis == 1) ? graph(smPt, 3) : graph(smPt, 2);
	  const auto l0i = l0*ndpc;
	  const auto r0i = r0*ndpc;

	  for (int i=0; i<ndpc; ++i){
	    m_memb.m_uMinusHalfNeg(i) = m_memb.m_state(l0i+i);
	    m_memb.m_uMinusHalfPos(i) = m_memb.m_state(uIndex+i);
	    m_memb.m_uPlusHalfNeg(i)  = m_memb.m_state(uIndex+i);
	    m_memb.m_uPlusHalfPos(i)  = m_memb.m_state(r0i+i);
	  }
	  break;
	}

      case InviscidFluxReconstruction::Weno3:
	{
	  const auto l0  = (axis == 1) ? graph(smPt, 1)  : graph(smPt, 4);
	  const auto r0  = (axis == 1) ? graph(smPt, 3)  : graph(smPt, 2);
	  const auto l1  = (axis == 1) ? graph(smPt, 5)  : graph(smPt, 8);
	  const auto r1  = (axis == 1) ? graph(smPt, 7)  : graph(smPt, 6);
	  const auto l0i = l0*ndpc;
	  const auto r0i = r0*ndpc;
	  const auto l1i = l1*ndpc;
	  const auto r1i = r1*ndpc;

	  pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg,
				 m_memb.m_uMinusHalfPos,
				 m_memb.m_uPlusHalfNeg,
				 m_memb.m_uPlusHalfPos,
				 m_memb.m_state,
				 l1i, l0i, r0i, r1i,
				 uIndex,
				 ndpc);
	  break;
	}

      case InviscidFluxReconstruction::Weno5:
	{
	  const auto l0  = (axis == 1) ? graph(smPt, 1)  : graph(smPt, 4);
	  const auto r0  = (axis == 1) ? graph(smPt, 3)  : graph(smPt, 2);
	  const auto l1  = (axis == 1) ? graph(smPt, 5)  : graph(smPt, 8);
	  const auto r1  = (axis == 1) ? graph(smPt, 7)  : graph(smPt, 6);
	  const auto l2  = (axis == 1) ? graph(smPt, 9)  : graph(smPt, 12);
	  const auto r2  = (axis == 1) ? graph(smPt, 11) : graph(smPt, 10);
	  const auto l0i = l0*ndpc;
	  const auto r0i = r0*ndpc;
	  const auto l1i = l1*ndpc;
	  const auto r1i = r1*ndpc;
	  const auto l2i = l2*ndpc;
	  const auto r2i = r2*ndpc;

	  pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg,
				 m_memb.m_uMinusHalfPos,
				 m_memb.m_uPlusHalfNeg,
				 m_memb.m_uPlusHalfPos,
				 m_memb.m_state,
				 l2i, l1i, l0i,
				 r0i, r1i, r2i,
				 uIndex,
				 ndpc);
	  break;
	}
      }
  }

private:
  using members_t = _ReconstructorMembers<edge_rec_t, state_t, mesh_t>;
  members_t m_memb;
};

// specialize for 3d
template<class edge_rec_t, class state_t, class mesh_t>
class ReconstructorFromState<3, edge_rec_t, state_t, mesh_t>
{

public:
  ReconstructorFromState() = delete;

  template<typename ...Args>
  ReconstructorFromState(Args && ... args)
    : m_memb(std::forward<Args>(args)...)
  {}

  template<int ndpc>
  void operator()(int smPt)
  {
    const auto & graph  = m_memb.m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0)*ndpc;
    const auto axis = m_memb.m_axis;

    /*
      For 3d, graph is:
      ptID l0 f0 r0 ba0 bot0 top0 [l1 f1 r1 ba1 bot1 top1 l2 f2 r2 ba2 bot2 top2]
       0    1  2  3  4    5   6   [ 7  8  9  10  11   12  13 14 15  16  17  18  ]

      where [] are optionally present, depending on the stencil chosen.
    */

    switch(m_memb.m_recEn)
      {
      case InviscidFluxReconstruction::FirstOrder:
	{
	  auto l0 = (axis == 1) ? graph(smPt, 1) : (axis==2) ? graph(smPt, 4) : graph(smPt, 5);
	  auto r0 = (axis == 1) ? graph(smPt, 3) : (axis==2) ? graph(smPt, 2) : graph(smPt, 6);
	  const auto l0i = l0*ndpc;
	  const auto r0i = r0*ndpc;

	  for (int i=0; i<ndpc; ++i){
	    m_memb.m_uMinusHalfNeg(i) = m_memb.m_state(l0i+i);
	    m_memb.m_uMinusHalfPos(i) = m_memb.m_state(uIndex+i);
	    m_memb.m_uPlusHalfNeg(i)  = m_memb.m_state(uIndex+i);
	    m_memb.m_uPlusHalfPos(i)  = m_memb.m_state(r0i+i);
	  }
	  break;
	}

      case InviscidFluxReconstruction::Weno3:
	{
	  auto l0 = (axis == 1) ? graph(smPt, 1) : (axis==2) ? graph(smPt, 4) : graph(smPt, 5);
	  auto r0 = (axis == 1) ? graph(smPt, 3) : (axis==2) ? graph(smPt, 2) : graph(smPt, 6);
	  auto l1 = (axis == 1) ? graph(smPt, 7) : (axis==2) ? graph(smPt, 10) : graph(smPt, 11);
	  auto r1 = (axis == 1) ? graph(smPt, 9) : (axis==2) ? graph(smPt, 8) : graph(smPt, 12);

	  const auto l0i = l0*ndpc;
	  const auto r0i = r0*ndpc;
	  const auto l1i = l1*ndpc;
	  const auto r1i = r1*ndpc;
	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg,
			       m_memb.m_uMinusHalfPos,
			       m_memb.m_uPlusHalfNeg,
			       m_memb.m_uPlusHalfPos,
			       m_memb.m_state,
			       l1i, l0i, r0i, r1i,
			       uIndex,
			       ndpc);
	  break;
	}

      case InviscidFluxReconstruction::Weno5:
	{
	  throw std::runtime_error("missing impl");
	  break;
	}
      }

  }

private:
  using members_t = _ReconstructorMembers<edge_rec_t, state_t, mesh_t>;
  members_t m_memb;
};

}}
#endif
