
#ifndef PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STATE_HPP_
#define PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STATE_HPP_

namespace pressiodemoapps{ namespace impl{

template<class edge_rec_t, class state_t, class mesh_t>
class ReconstructorFromState
{

public:
  ReconstructorFromState() = delete;
  ReconstructorFromState(const int axis,
			 ::pressiodemoapps::reconstructionEnum recEn,
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

  template<int ndpc>
  void operator()(int smPt)
  {
    const auto & graph = m_meshObj.graph();
    const auto uIndex = graph(smPt, 0)*ndpc;

    switch(m_recEn){
    case pressiodemoapps::reconstructionEnum::firstOrder:
      {
	const auto l0 = (m_axis == 1) ? graph(smPt, 1) : graph(smPt, 4);
	const auto r0 = (m_axis == 1) ? graph(smPt, 3) : graph(smPt, 2);
	const auto l0i = l0*ndpc;
	const auto r0i = r0*ndpc;

	for (int i=0; i<ndpc; ++i){
	  m_uMinusHalfNeg(i) = m_state(l0i+i);
	  m_uMinusHalfPos(i) = m_state(uIndex+i);
	  m_uPlusHalfNeg(i)  = m_state(uIndex+i);
	  m_uPlusHalfPos(i)  = m_state(r0i+i);
	}
	break;
      }

    case pressiodemoapps::reconstructionEnum::fifthOrderWeno:
      {
	const auto l0 = (m_axis == 1) ? graph(smPt, 1)  : graph(smPt, 4);
	const auto r0 = (m_axis == 1) ? graph(smPt, 3)  : graph(smPt, 2);
	const auto l1 = (m_axis == 1) ? graph(smPt, 5)  : graph(smPt, 8);
	const auto r1 = (m_axis == 1) ? graph(smPt, 7)  : graph(smPt, 6);
	const auto l2 = (m_axis == 1) ? graph(smPt, 9)  : graph(smPt, 12);
	const auto r2 = (m_axis == 1) ? graph(smPt, 11) : graph(smPt, 10);
	const auto l0i = l0*ndpc;
	const auto r0i = r0*ndpc;
	const auto l1i = l1*ndpc;
	const auto r1i = r1*ndpc;
	const auto l2i = l2*ndpc;
	const auto r2i = r2*ndpc;

	pressiodemoapps::weno5(m_uMinusHalfNeg,
			       m_uMinusHalfPos,
			       m_uPlusHalfNeg,
			       m_uPlusHalfPos,
			       m_state,
			       l2i, l1i, l0i,
			       r0i, r1i, r2i,
			       uIndex,
			       ndpc);
	break;
      }
    }
  }

private:
  int m_axis;
  ::pressiodemoapps::reconstructionEnum m_recEn;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  edge_rec_t & m_uMinusHalfNeg;
  edge_rec_t & m_uMinusHalfPos;
  edge_rec_t & m_uPlusHalfNeg;
  edge_rec_t & m_uPlusHalfPos;

};

}}
#endif
