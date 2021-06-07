
#ifndef PRESSIODEMOAPPS_GHOST_FILLER_EE2D_HPP_
#define PRESSIODEMOAPPS_GHOST_FILLER_EE2D_HPP_

namespace pressiodemoapps{ namespace impl{

template<
  int numDofPerCell,
  class state_t,
  class mesh_t,
  class ghost_t
  >
class Ghost2dNeumannFiller;

template<class state_t, class mesh_t, class ghost_t>
class Ghost2dNeumannFiller<4, state_t, mesh_t, ghost_t>
{

public:
  Ghost2dNeumannFiller() = delete;
  Ghost2dNeumannFiller(const int stencilSize,
		       const state_t & stateIn,
		       const mesh_t & meshIn,
		       ghost_t & ghostLeft,
		       ghost_t & ghostTop,
		       ghost_t & ghostRight,
		       ghost_t & ghostBottom)
    : m_stencilSize(stencilSize),
      m_state(stateIn),
      m_meshObj(meshIn),
      m_ghostLeft(ghostLeft),
      m_ghostTop(ghostTop),
      m_ghostRight(ghostRight),
      m_ghostBottom(ghostBottom)
  {}

  template<class index_t>
  void operator()(index_t smPt)
  {
    const auto & graph = m_meshObj.graph();
    const auto cellGID = graph(smPt, 0);
    const auto w0 = graph(smPt, 1);
    const auto n0 = graph(smPt, 2);
    const auto e0 = graph(smPt, 3);
    const auto s0 = graph(smPt, 4);

    if (w0 == -1)
    {
      const auto ind = cellGID*4;
      m_ghostLeft(smPt, 0) = m_state(ind);
      m_ghostLeft(smPt, 1) = m_state(ind+1);
      m_ghostLeft(smPt, 2) = m_state(ind+2);
      m_ghostLeft(smPt, 3) = m_state(ind+3);
    }

    if (n0 == -1)
    {
      const auto ind = cellGID*4;
      m_ghostTop(smPt, 0) = m_state(ind);
      m_ghostTop(smPt, 1) = m_state(ind+1);
      m_ghostTop(smPt, 2) = m_state(ind+2);
      m_ghostTop(smPt, 3) = m_state(ind+3);
    }

    if (e0 == -1)
    {
      const auto ind = cellGID*4;
      m_ghostRight(smPt, 0) = m_state(ind);
      m_ghostRight(smPt, 1) = m_state(ind+1);
      m_ghostRight(smPt, 2) = m_state(ind+2);
      m_ghostRight(smPt, 3) = m_state(ind+3);
    }

    if (s0 == -1)
    {
      const auto ind = cellGID*4;
      m_ghostBottom(smPt, 0) = m_state(ind);
      m_ghostBottom(smPt, 1) = m_state(ind+1);
      m_ghostBottom(smPt, 2) = m_state(ind+2);
      m_ghostBottom(smPt, 3) = m_state(ind+3);
    }

    if (m_stencilSize == 7)
      {
	const auto w1 = graph(smPt, 5);
	const auto n1 = graph(smPt, 6);
	const auto e1 = graph(smPt, 7);
	const auto s1 = graph(smPt, 8);
	const auto w2 = graph(smPt, 9);
	const auto n2 = graph(smPt, 10);
	const auto e2 = graph(smPt, 11);
	const auto s2 = graph(smPt, 12);

	if (w1 == -1){
	  const auto ind = e0*4;
	  m_ghostLeft(smPt, 4) = m_state(ind);
	  m_ghostLeft(smPt, 5) = m_state(ind+1);
	  m_ghostLeft(smPt, 6) = m_state(ind+2);
	  m_ghostLeft(smPt, 7) = m_state(ind+3);
	}

	if (n1 == -1){
	  const auto ind = s0*4;
	  m_ghostTop(smPt, 4) = m_state(ind);
	  m_ghostTop(smPt, 5) = m_state(ind+1);
	  m_ghostTop(smPt, 6) = m_state(ind+2);
	  m_ghostTop(smPt, 7) = m_state(ind+3);
	}

	if (e1 == -1){
	  const auto ind = w0*4;
	  m_ghostRight(smPt, 4) = m_state(ind);
	  m_ghostRight(smPt, 5) = m_state(ind+1);
	  m_ghostRight(smPt, 6) = m_state(ind+2);
	  m_ghostRight(smPt, 7) = m_state(ind+3);
	}

	if (s1 == -1){
	  const auto ind = n0*4;
	  m_ghostBottom(smPt, 4) = m_state(ind);
	  m_ghostBottom(smPt, 5) = m_state(ind+1);
	  m_ghostBottom(smPt, 6) = m_state(ind+2);
	  m_ghostBottom(smPt, 7) = m_state(ind+3);
	}

	if (w2 == -1){
	  const auto ind = e1*4;
	  m_ghostLeft(smPt, 8)  = m_state(ind);
	  m_ghostLeft(smPt, 9)  = m_state(ind+1);
	  m_ghostLeft(smPt, 10) = m_state(ind+2);
	  m_ghostLeft(smPt, 11) = m_state(ind+3);
	}

	if (n2 == -1){
	  const auto ind = s1*4;
	  m_ghostTop(smPt, 8) = m_state(ind);
	  m_ghostTop(smPt, 9) = m_state(ind+1);
	  m_ghostTop(smPt, 10) = m_state(ind+2);
	  m_ghostTop(smPt, 11) = m_state(ind+3);
	}

	if (e2 == -1){
	  const auto ind = w1*4;
	  m_ghostRight(smPt, 8) = m_state(ind);
	  m_ghostRight(smPt, 9) = m_state(ind+1);
	  m_ghostRight(smPt, 10) = m_state(ind+2);
	  m_ghostRight(smPt, 11) = m_state(ind+3);
	}

	if (s2 == -1){
	  const auto ind = n1*4;
	  m_ghostBottom(smPt, 8) = m_state(ind);
	  m_ghostBottom(smPt, 9) = m_state(ind+1);
	  m_ghostBottom(smPt, 10) = m_state(ind+2);
	  m_ghostBottom(smPt, 11) = m_state(ind+3);
	}
       }

  }

private:
  const int m_stencilSize;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  ghost_t & m_ghostLeft;
  ghost_t & m_ghostTop;
  ghost_t & m_ghostRight;
  ghost_t & m_ghostBottom;
};

}}
#endif
