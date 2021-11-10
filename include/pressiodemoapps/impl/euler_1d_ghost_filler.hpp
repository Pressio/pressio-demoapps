
#ifndef PRESSIODEMOAPPS_GHOST_FILLER_EE1D_HPP_
#define PRESSIODEMOAPPS_GHOST_FILLER_EE1D_HPP_

namespace pressiodemoapps{ namespace impl{

template<
  int numDofPerCell,
  class state_t,
  class mesh_t,
  class ghost_t
  >
class Ghost1dNeumannFiller;

template<class state_t, class mesh_t, class ghost_t>
class Ghost1dNeumannFiller<3, state_t, mesh_t, ghost_t>
{

public:
  Ghost1dNeumannFiller() = delete;
  Ghost1dNeumannFiller(const int stencilSize,
		       const state_t & stateIn,
		       const mesh_t & meshIn,
		       ghost_t & ghostLeft,
		       ghost_t & ghostRight)
    : m_stencilSize(stencilSize),
      m_state(stateIn),
      m_meshObj(meshIn),
      m_ghostLeft(ghostLeft),
      m_ghostRight(ghostRight)
  {}

  template<class index_t>
  void operator()(index_t smPt, int gRow)
  {
    constexpr int numDofPerCell = 3;
    const auto & graph = m_meshObj.graph();
    const auto cellGID = graph(smPt, 0);
    const auto uIndex  = cellGID*numDofPerCell;

    const auto left0  = graph(smPt, 1);
    const auto right0 = graph(smPt, 2);

    if (left0 == -1){
      m_ghostLeft(gRow, 0) = m_state(uIndex);
      m_ghostLeft(gRow, 1) = m_state(uIndex+1);
      m_ghostLeft(gRow, 2) = m_state(uIndex+2);
    }

    if (right0 == -1){
      m_ghostRight(gRow, 0) = m_state(uIndex);
      m_ghostRight(gRow, 1) = m_state(uIndex+1);
      m_ghostRight(gRow, 2) = m_state(uIndex+2);
    }

    if (m_stencilSize>=5)
      {
	const auto left1  = graph(smPt, 3);
	const auto right1 = graph(smPt, 4);

	if (left1 == -1){
	  const auto ind = right0*numDofPerCell;
	  m_ghostLeft(gRow, 3) = m_state(ind);
	  m_ghostLeft(gRow, 4) = m_state(ind+1);
	  m_ghostLeft(gRow, 5) = m_state(ind+2);
	}

	if (right1 == -1){
	  const auto ind = left0*numDofPerCell;
	  m_ghostRight(gRow, 3) = m_state(ind);
	  m_ghostRight(gRow, 4) = m_state(ind+1);
	  m_ghostRight(gRow, 5) = m_state(ind+2);
	}
      }

    if (m_stencilSize == 7){
      const auto left1  = graph(smPt, 5);
      const auto right1 = graph(smPt, 6);
      const auto left2  = graph(smPt, 7);
      const auto right2 = graph(smPt, 8);

      if (left2 == -1){
	const auto ind = right1*numDofPerCell;
	m_ghostLeft(gRow, 6) = m_state(ind);
	m_ghostLeft(gRow, 7) = m_state(ind+1);
	m_ghostLeft(gRow, 8) = m_state(ind+2);
      }

      if (right2 == -1){
	const auto ind = left1*numDofPerCell;
	m_ghostRight(gRow, 6) = m_state(ind);
	m_ghostRight(gRow, 7) = m_state(ind+1);
	m_ghostRight(gRow, 8) = m_state(ind+2);
      }
    }
  }

private:
  const int m_stencilSize;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  ghost_t & m_ghostLeft;
  ghost_t & m_ghostRight;
};

}}
#endif
