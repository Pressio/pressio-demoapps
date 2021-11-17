
#ifndef PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STATE_HPP_
#define PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STATE_HPP_

namespace pressiodemoapps{ namespace impl{

template<class ...Args>
struct _ReconstructorMembers;

// specialize for VOID reconstruction gradients
template<class edge_rec_t, class state_t, class mesh_t>
struct _ReconstructorMembers<edge_rec_t, state_t, mesh_t>
{
  _ReconstructorMembers() = delete;

  _ReconstructorMembers(const int axis, InviscidFluxReconstruction recEn,
			const state_t & state, const mesh_t & meshObj,
			edge_rec_t & uMinusHalfNeg, edge_rec_t & uMinusHalfPos,
			edge_rec_t & uPlusHalfNeg,  edge_rec_t & uPlusHalfPos)
    : m_axis(axis), m_recEn(recEn), m_state(state), m_meshObj(meshObj),
      m_uMinusHalfNeg(uMinusHalfNeg), m_uMinusHalfPos(uMinusHalfPos),
      m_uPlusHalfNeg(uPlusHalfNeg),   m_uPlusHalfPos(uPlusHalfPos)
  {}

  int m_axis = -1;
  InviscidFluxReconstruction m_recEn;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  edge_rec_t & m_uMinusHalfNeg;
  edge_rec_t & m_uMinusHalfPos;
  edge_rec_t & m_uPlusHalfNeg;
  edge_rec_t & m_uPlusHalfPos;
};

// reconstruction gradients are needed
template<class edge_rec_t, class state_t, class mesh_t, class rec_grad_t>
struct _ReconstructorMembers<edge_rec_t, state_t, mesh_t, rec_grad_t>
{
  _ReconstructorMembers() = delete;

  _ReconstructorMembers(const int axis, InviscidFluxReconstruction recEn,
			const state_t & state, const mesh_t & meshObj,
			edge_rec_t & uMinusHalfNeg, edge_rec_t & uMinusHalfPos,
			edge_rec_t & uPlusHalfNeg,  edge_rec_t & uPlusHalfPos,
			rec_grad_t & gradLNeg, rec_grad_t & gradLPos,
			rec_grad_t & gradRNeg, rec_grad_t & gradRPos)
    : m_recEn(recEn), m_state(state), m_meshObj(meshObj),
      m_uMinusHalfNeg(uMinusHalfNeg), m_uMinusHalfPos(uMinusHalfPos),
      m_uPlusHalfNeg(uPlusHalfNeg),   m_uPlusHalfPos(uPlusHalfPos),
      m_gradLNeg(gradLNeg), m_gradLPos(gradLPos),
      m_gradRNeg(gradRNeg), m_gradRPos(gradRPos)
  {}

  int m_axis = -1;
  InviscidFluxReconstruction m_recEn;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  edge_rec_t & m_uMinusHalfNeg;
  edge_rec_t & m_uMinusHalfPos;
  edge_rec_t & m_uPlusHalfNeg;
  edge_rec_t & m_uPlusHalfPos;
  rec_grad_t & m_gradLNeg ;
  rec_grad_t & m_gradLPos ;
  rec_grad_t & m_gradRNeg ;
  rec_grad_t & m_gradRPos ;
};



// class edge_rec_t, class state_t, class mesh_t, class rec_grad_t >
template<int dim, int numDofPerCell, class ...Args>
class ReconstructorFromState;

// -----------------------------------------------------------
// 1d, multiple dofsPerCell, no reconstruction gradients
// -----------------------------------------------------------
template<int ndpc, class edge_rec_t, class state_t, class mesh_t>
class ReconstructorFromState<1, ndpc, edge_rec_t, state_t, mesh_t>
{
public:
  template<typename ...Args>
  ReconstructorFromState(Args && ... args) : m_memb(std::forward<Args>(args)...){}

  template<class IndexType>
  void operator()(IndexType smPt)
  {
    const auto & graph  = m_memb.m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0)*ndpc;
    const auto axis = m_memb.m_axis;

    /* For 1d, graph is:
       ptID left0 right0 [left1 right1 left2 right2]
        0     1     2    [  3     4      5     6   ]

      where [] are optionally present, depending on the stencil chosen.
    */

    const auto l0i = graph(smPt, 1)*ndpc;
    const auto r0i = graph(smPt, 2)*ndpc;

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
	  const auto l1i  = graph(smPt, 3)*ndpc;
	  const auto r1i  = graph(smPt, 4)*ndpc;
	  ::pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg,
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
	  const auto l1i = graph(smPt, 3)*ndpc;
	  const auto r1i = graph(smPt, 4)*ndpc;
	  const auto l2i = graph(smPt, 5)*ndpc;
	  const auto r2i = graph(smPt, 6)*ndpc;
	  ::pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg,
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

// -----------------------------------------------------------
// 1d, multiple dofsPerCell, WITH reconstruction gradients
// -----------------------------------------------------------
template<int ndpc, class edge_rec_t, class state_t, class mesh_t, class rec_grad_t>
class ReconstructorFromState<1, ndpc, edge_rec_t, state_t, mesh_t, rec_grad_t>
{
public:
  template<typename ...Args>
  ReconstructorFromState(Args && ... args) : m_memb(std::forward<Args>(args)...){}

  template<class IndexType>
  void operator()(IndexType smPt)
  {
    const auto & graph  = m_memb.m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0)*ndpc;
    const auto axis = m_memb.m_axis;

    const auto l0i = graph(smPt, 1)*ndpc;
    const auto r0i = graph(smPt, 2)*ndpc;

    switch(m_memb.m_recEn)
      {
      case InviscidFluxReconstruction::FirstOrder:
	{
	  for (int i=0; i<ndpc; ++i){
	    m_memb.m_uMinusHalfNeg(i) = m_memb.m_state(l0i+i);
	    m_memb.m_uMinusHalfPos(i) = m_memb.m_state(uIndex+i);
	    m_memb.m_uPlusHalfNeg(i)  = m_memb.m_state(uIndex+i);
	    m_memb.m_uPlusHalfPos(i)  = m_memb.m_state(r0i+i);
	    m_memb.m_gradLNeg(i, 0) = 1;
	    m_memb.m_gradLNeg(i, 1) = 0;
	    m_memb.m_gradLPos(i, 0) = 0;
	    m_memb.m_gradLPos(i, 1) = 1;
	    m_memb.m_gradRNeg(i, 0) = 1;
	    m_memb.m_gradRNeg(i, 1) = 0;
	    m_memb.m_gradRPos(i, 0) = 0;
	    m_memb.m_gradRPos(i, 1) = 1;
	  }

	  break;
	}

      case InviscidFluxReconstruction::Weno3:
	{
	  const auto l1i  = graph(smPt, 3)*ndpc;
	  const auto r1i  = graph(smPt, 4)*ndpc;
	  ::pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg,
				   m_memb.m_uMinusHalfPos,
				   m_memb.m_uPlusHalfNeg,
				   m_memb.m_uPlusHalfPos,
				   m_memb.m_gradLNeg,
				   m_memb.m_gradLPos,
				   m_memb.m_gradRNeg,
				   m_memb.m_gradRPos,
				   m_memb.m_state,
				   l1i, l0i, r0i, r1i,
				   uIndex,
				   ndpc);
	  break;
	}

      case InviscidFluxReconstruction::Weno5:
	{
	  const auto l1i = graph(smPt, 3)*ndpc;
	  const auto r1i = graph(smPt, 4)*ndpc;
	  const auto l2i = graph(smPt, 5)*ndpc;
	  const auto r2i = graph(smPt, 6)*ndpc;
	  ::pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg,
				   m_memb.m_uMinusHalfPos,
				   m_memb.m_uPlusHalfNeg,
				   m_memb.m_uPlusHalfPos,
				   m_memb.m_gradLNeg,
				   m_memb.m_gradLPos,
				   m_memb.m_gradRNeg,
				   m_memb.m_gradRPos,
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
  using members_t = _ReconstructorMembers<
  edge_rec_t, state_t, mesh_t, rec_grad_t>;
  members_t m_memb;
};


// -----------------------------------------------------------
// 1d, 1 dof, no reconstruction gradients
// -----------------------------------------------------------
template<class edge_rec_t, class state_t, class mesh_t>
class ReconstructorFromState<1, 1, edge_rec_t, state_t, mesh_t>{
public:
  template<typename ...Args>
  ReconstructorFromState(Args && ... args) : m_memb(std::forward<Args>(args)...){}

  template<class IndexType>
  void operator()(IndexType smPt)
  {
    const auto & graph  = m_memb.m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0);
    const auto axis = m_memb.m_axis;

    const auto l0i = graph(smPt, 1);
    const auto r0i = graph(smPt, 2);
    switch(m_memb.m_recEn)
      {
      case InviscidFluxReconstruction::FirstOrder:
	{
	  m_memb.m_uMinusHalfNeg = m_memb.m_state(l0i);
	  m_memb.m_uMinusHalfPos = m_memb.m_state(uIndex);
	  m_memb.m_uPlusHalfNeg  = m_memb.m_state(uIndex);
	  m_memb.m_uPlusHalfPos  = m_memb.m_state(r0i);
	  break;
	}

      case InviscidFluxReconstruction::Weno3:
	{
	  const auto l1i = graph(smPt, 3);
	  const auto r1i = graph(smPt, 4);
	  ::pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg,
				   m_memb.m_uMinusHalfPos,
				   m_memb.m_uPlusHalfNeg,
				   m_memb.m_uPlusHalfPos,
				   m_memb.m_state(l1i),
				   m_memb.m_state(l0i),
				   m_memb.m_state(uIndex),
				   m_memb.m_state(r0i),
				   m_memb.m_state(r1i));
	  break;
	}

      case InviscidFluxReconstruction::Weno5:
	{
	  const auto l1i = graph(smPt, 3);
	  const auto r1i = graph(smPt, 4);
	  const auto l2i = graph(smPt, 5);
	  const auto r2i = graph(smPt, 6);
	  ::pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg,
				   m_memb.m_uMinusHalfPos,
				   m_memb.m_uPlusHalfNeg,
				   m_memb.m_uPlusHalfPos,
				   m_memb.m_state(l2i),
				   m_memb.m_state(l1i),
				   m_memb.m_state(l0i),
				   m_memb.m_state(uIndex),
				   m_memb.m_state(r0i),
				   m_memb.m_state(r1i),
				   m_memb.m_state(r2i));
	  break;
	}
      }
  }

private:
  using members_t = _ReconstructorMembers<edge_rec_t, state_t, mesh_t>;
  members_t m_memb;
};

// -----------------------------------------------------------
// 1d, 1 dofsPerCell, WITH reconstruction gradients
// -----------------------------------------------------------
template<class edge_rec_t, class state_t, class mesh_t, class rec_grad_t>
class ReconstructorFromState<1, 1, edge_rec_t, state_t, mesh_t, rec_grad_t>{
public:
  template<typename ...Args>
  ReconstructorFromState(Args && ... args) : m_memb(std::forward<Args>(args)...){}

  template<class IndexType>
  void operator()(IndexType smPt)
  {
    const auto & graph  = m_memb.m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0);
    const auto axis = m_memb.m_axis;

    const auto l0i  = graph(smPt, 1);
    const auto r0i  = graph(smPt, 2);
    switch(m_memb.m_recEn)
      {
      case InviscidFluxReconstruction::FirstOrder:
	{
	  m_memb.m_uMinusHalfNeg = m_memb.m_state(l0i);
	  m_memb.m_uMinusHalfPos = m_memb.m_state(uIndex);
	  m_memb.m_uPlusHalfNeg  = m_memb.m_state(uIndex);
	  m_memb.m_uPlusHalfPos  = m_memb.m_state(r0i);
	  m_memb.m_gradLNeg[0] = 1;
	  m_memb.m_gradLNeg[1] = 0;
	  m_memb.m_gradLPos[0] = 0;
	  m_memb.m_gradLPos[1] = 1;
	  m_memb.m_gradRNeg[0] = 1;
	  m_memb.m_gradRNeg[1] = 0;
	  m_memb.m_gradRPos[0] = 0;
	  m_memb.m_gradRPos[1] = 1;
	  break;
	}

      case InviscidFluxReconstruction::Weno3:
	{
	  const auto l1i = graph(smPt, 3);
	  const auto r1i = graph(smPt, 4);
	  ::pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg,
				   m_memb.m_uMinusHalfPos,
				   m_memb.m_uPlusHalfNeg,
				   m_memb.m_uPlusHalfPos,
				   m_memb.m_gradLNeg,
				   m_memb.m_gradLPos,
				   m_memb.m_gradRNeg,
				   m_memb.m_gradRPos,
				   m_memb.m_state(l1i),
				   m_memb.m_state(l0i),
				   m_memb.m_state(uIndex),
				   m_memb.m_state(r0i),
				   m_memb.m_state(r1i));
	  break;
	}

      case InviscidFluxReconstruction::Weno5:
	{
	  const auto l1i  = graph(smPt, 3);
	  const auto r1i  = graph(smPt, 4);
	  const auto l2i  = graph(smPt, 5);
	  const auto r2i  = graph(smPt, 6);
	  ::pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg,
				   m_memb.m_uMinusHalfPos,
				   m_memb.m_uPlusHalfNeg,
				   m_memb.m_uPlusHalfPos,
				   m_memb.m_gradLNeg,
				   m_memb.m_gradLPos,
				   m_memb.m_gradRNeg,
				   m_memb.m_gradRPos,
				   m_memb.m_state(l2i),
				   m_memb.m_state(l1i),
				   m_memb.m_state(l0i),
				   m_memb.m_state(uIndex),
				   m_memb.m_state(r0i),
				   m_memb.m_state(r1i),
				   m_memb.m_state(r2i));
	  break;
	}
      }
  }

private:
  using members_t = _ReconstructorMembers<edge_rec_t, state_t, mesh_t, rec_grad_t>;
  members_t m_memb;
};


// -----------------------------------------------------------
// 2d, multiple dofsPerCell, no reconstruction gradients
// -----------------------------------------------------------
template<int ndpc, class edge_rec_t, class state_t, class mesh_t>
class ReconstructorFromState<2, ndpc, edge_rec_t, state_t, mesh_t>{
public:
  template<typename ...Args>
  ReconstructorFromState(Args && ... args) : m_memb(std::forward<Args>(args)...){}

  template<class IndexType>
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

	  ::pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg,
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

	  ::pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg,
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


// -----------------------------------------------------------
// 2d, multiple dofsPerCell, WITH reconstruction gradients
// -----------------------------------------------------------
template<int ndpc, class edge_rec_t, class state_t, class mesh_t, class rec_grad_t>
class ReconstructorFromState<2, ndpc, edge_rec_t, state_t, mesh_t, rec_grad_t>{
public:
  template<typename ...Args>
  ReconstructorFromState(Args && ... args) : m_memb(std::forward<Args>(args)...){}

  template<class IndexType>
  void operator()(IndexType smPt)
  {
    const auto & graph  = m_memb.m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0)*ndpc;
    const auto axis = m_memb.m_axis;

    switch(m_memb.m_recEn)
      {
      case InviscidFluxReconstruction::FirstOrder:
	{
	  throw std::runtime_error("REC MISSING C");
	  // const auto l0  = (axis == 1) ? graph(smPt, 1) : graph(smPt, 4);
	  // const auto r0  = (axis == 1) ? graph(smPt, 3) : graph(smPt, 2);
	  // const auto l0i = l0*ndpc;
	  // const auto r0i = r0*ndpc;

	  // for (int i=0; i<ndpc; ++i){
	  //   m_memb.m_uMinusHalfNeg(i) = m_memb.m_state(l0i+i);
	  //   m_memb.m_uMinusHalfPos(i) = m_memb.m_state(uIndex+i);
	  //   m_memb.m_uPlusHalfNeg(i)  = m_memb.m_state(uIndex+i);
	  //   m_memb.m_uPlusHalfPos(i)  = m_memb.m_state(r0i+i);
	  // }
	  break;
	}

      case InviscidFluxReconstruction::Weno3:
	{
	  throw std::runtime_error("REC MISSING D");
	  // const auto l0  = (axis == 1) ? graph(smPt, 1)  : graph(smPt, 4);
	  // const auto r0  = (axis == 1) ? graph(smPt, 3)  : graph(smPt, 2);
	  // const auto l1  = (axis == 1) ? graph(smPt, 5)  : graph(smPt, 8);
	  // const auto r1  = (axis == 1) ? graph(smPt, 7)  : graph(smPt, 6);
	  // const auto l0i = l0*ndpc;
	  // const auto r0i = r0*ndpc;
	  // const auto l1i = l1*ndpc;
	  // const auto r1i = r1*ndpc;

	  // ::pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg,
	  // 			 m_memb.m_uMinusHalfPos,
	  // 			 m_memb.m_uPlusHalfNeg,
	  // 			 m_memb.m_uPlusHalfPos,
	  // 			 m_memb.m_state,
	  // 			 l1i, l0i, r0i, r1i,
	  // 			 uIndex,
	  // 			 ndpc);
	  break;
	}

      case InviscidFluxReconstruction::Weno5:
	{
	  throw std::runtime_error("REC MISSING E");
	  // const auto l0  = (axis == 1) ? graph(smPt, 1)  : graph(smPt, 4);
	  // const auto r0  = (axis == 1) ? graph(smPt, 3)  : graph(smPt, 2);
	  // const auto l1  = (axis == 1) ? graph(smPt, 5)  : graph(smPt, 8);
	  // const auto r1  = (axis == 1) ? graph(smPt, 7)  : graph(smPt, 6);
	  // const auto l2  = (axis == 1) ? graph(smPt, 9)  : graph(smPt, 12);
	  // const auto r2  = (axis == 1) ? graph(smPt, 11) : graph(smPt, 10);
	  // const auto l0i = l0*ndpc;
	  // const auto r0i = r0*ndpc;
	  // const auto l1i = l1*ndpc;
	  // const auto r1i = r1*ndpc;
	  // const auto l2i = l2*ndpc;
	  // const auto r2i = r2*ndpc;

	  // ::pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg,
	  // 			 m_memb.m_uMinusHalfPos,
	  // 			 m_memb.m_uPlusHalfNeg,
	  // 			 m_memb.m_uPlusHalfPos,
	  // 			 m_memb.m_state,
	  // 			 l2i, l1i, l0i,
	  // 			 r0i, r1i, r2i,
	  // 			 uIndex,
	  // 			 ndpc);
	  break;
	}
      }
  }

private:
  using members_t = _ReconstructorMembers<edge_rec_t, state_t, mesh_t, rec_grad_t>;
  members_t m_memb;
};


// -----------------------------------------------------------
// 3d, multiple dofsPerCell, no reconstruction gradients
// -----------------------------------------------------------
template<int ndpc, class edge_rec_t, class state_t, class mesh_t>
class ReconstructorFromState<3, ndpc, edge_rec_t, state_t, mesh_t>{

public:
  template<typename ...Args>
  ReconstructorFromState(Args && ... args) : m_memb(std::forward<Args>(args)...){}

  template<class IndexType>
  void operator()(IndexType smPt)
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
	::pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg,
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
