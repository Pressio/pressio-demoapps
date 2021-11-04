
#ifndef PRESSIODEMOAPPS_GHOST_FILLER_EE2D_SEDOV2D_SYM_HPP_
#define PRESSIODEMOAPPS_GHOST_FILLER_EE2D_SEDOV2D_SYM_HPP_

namespace pressiodemoapps{ namespace impl{

template<class state_t, class mesh_t, class ghost_t>
class Sedov2dSymmetryGhostFiller
{

public:
  Sedov2dSymmetryGhostFiller() = delete;
  Sedov2dSymmetryGhostFiller(const int stencilSize,
			     const state_t & stateIn,
			     const mesh_t & meshIn,
			     ghost_t & ghostLeft,
			     ghost_t & ghostFront,
			     ghost_t & ghostRight,
			     ghost_t & ghostBack)
    : m_stencilSize(stencilSize),
      m_state(stateIn),
      m_meshObj(meshIn),
      m_ghostLeft(ghostLeft),
      m_ghostFront(ghostFront),
      m_ghostRight(ghostRight),
      m_ghostBack(ghostBack)
  {}

  template<int stencilSize, class index_t>
  typename std::enable_if<stencilSize == 3>::type
  operator()(index_t smPt, int gRow)
  {

    constexpr int numDofPerCell = 4;
    const auto & graph = m_meshObj.graph();
    assert(::pressiodemoapps::extent(graph, 0) >= 5);
    const auto cellGID = graph(smPt, 0);
    const auto uIndex  = cellGID*numDofPerCell;

    const auto left0 = graph(smPt, 1);
    const auto front0 = graph(smPt, 2);
    const auto right0 = graph(smPt, 3);
    const auto back0 = graph(smPt, 4);

    if (left0 == -1)
    {
      m_ghostLeft(gRow, 0) = m_state(uIndex);
      m_ghostLeft(gRow, 1) = -m_state(uIndex+1);
      m_ghostLeft(gRow, 2) = m_state(uIndex+2);
      m_ghostLeft(gRow, 3) = m_state(uIndex+3);
    }

    if (front0 == -1)
    {
      m_ghostFront(gRow, 0) = m_state(uIndex);
      m_ghostFront(gRow, 1) = m_state(uIndex+1);
      m_ghostFront(gRow, 2) = m_state(uIndex+2);
      m_ghostFront(gRow, 3) = m_state(uIndex+3);
    }

    if (right0 == -1)
    {
      m_ghostRight(gRow, 0) = m_state(uIndex);
      m_ghostRight(gRow, 1) = m_state(uIndex+1);
      m_ghostRight(gRow, 2) = m_state(uIndex+2);
      m_ghostRight(gRow, 3) = m_state(uIndex+3);
    }

    if (back0 == -1)
    {
      m_ghostBack(gRow, 0) = m_state(uIndex);
      m_ghostBack(gRow, 1) = m_state(uIndex+1);
      m_ghostBack(gRow, 2) = -m_state(uIndex+2);
      m_ghostBack(gRow, 3) = m_state(uIndex+3);
    }
  }

  template<int stencilSize, class index_t>
  typename std::enable_if<stencilSize == 5>::type
  operator()(index_t smPt, int gRow)
  {
    this->template operator()<3, index_t>(smPt, gRow);

    constexpr int numDofPerCell = 4;
    const auto & graph = m_meshObj.graph();
    assert(::pressiodemoapps::extent(graph, 0) >= 9);

    const auto left0 = graph(smPt, 1);
    const auto front0 = graph(smPt, 2);
    const auto right0 = graph(smPt, 3);
    const auto back0 = graph(smPt, 4);
    const auto left1 = graph(smPt, 5);
    const auto front1 = graph(smPt, 6);
    const auto right1 = graph(smPt, 7);
    const auto back1 = graph(smPt, 8);

    if (left1 == -1){
      const auto ind = right0*numDofPerCell;
      m_ghostLeft(gRow, 4) = m_state(ind);
      m_ghostLeft(gRow, 5) = -m_state(ind+1);
      m_ghostLeft(gRow, 6) = m_state(ind+2);
      m_ghostLeft(gRow, 7) = m_state(ind+3);
    }

    if (front1 == -1){
      const auto ind = back0*numDofPerCell;
      m_ghostFront(gRow, 4) = m_state(ind);
      m_ghostFront(gRow, 5) = m_state(ind+1);
      m_ghostFront(gRow, 6) = m_state(ind+2);
      m_ghostFront(gRow, 7) = m_state(ind+3);
    }

    if (right1 == -1){
      const auto ind = left0*numDofPerCell;
      m_ghostRight(gRow, 4) = m_state(ind);
      m_ghostRight(gRow, 5) = m_state(ind+1);
      m_ghostRight(gRow, 6) = m_state(ind+2);
      m_ghostRight(gRow, 7) = m_state(ind+3);
    }

    if (back1 == -1){
      const auto ind = front0*numDofPerCell;
      m_ghostBack(gRow, 4) = m_state(ind);
      m_ghostBack(gRow, 5) = m_state(ind+1);
      m_ghostBack(gRow, 6) = -m_state(ind+2);
      m_ghostBack(gRow, 7) = m_state(ind+3);
    }
  }

  template<int stencilSize, class index_t>
  typename std::enable_if<stencilSize == 7>::type
  operator()(index_t smPt, int gRow)
  {
    this->template operator()<5, index_t>(smPt, gRow);

    constexpr int numDofPerCell = 4;
    const auto & graph = m_meshObj.graph();
    assert(::pressiodemoapps::extent(graph, 0) >= 13);

    const auto left1 = graph(smPt, 5);
    const auto front1 = graph(smPt, 6);
    const auto right1 = graph(smPt, 7);
    const auto back1 = graph(smPt, 8);
    const auto left2 = graph(smPt, 9);
    const auto front2 = graph(smPt, 10);
    const auto right2 = graph(smPt, 11);
    const auto back2 = graph(smPt, 12);

    if (left2 == -1){
      const auto ind = right1*numDofPerCell;
      m_ghostLeft(gRow, 8)  = m_state(ind);
      m_ghostLeft(gRow, 9)  = -m_state(ind+1);
      m_ghostLeft(gRow, 10) = m_state(ind+2);
      m_ghostLeft(gRow, 11) = m_state(ind+3);
    }

    if (front2 == -1){
      const auto ind = back1*numDofPerCell;
      m_ghostFront(gRow, 8) = m_state(ind);
      m_ghostFront(gRow, 9) = m_state(ind+1);
      m_ghostFront(gRow, 10) = m_state(ind+2);
      m_ghostFront(gRow, 11) = m_state(ind+3);
    }

    if (right2 == -1){
      const auto ind = left1*numDofPerCell;
      m_ghostRight(gRow, 8) = m_state(ind);
      m_ghostRight(gRow, 9) = m_state(ind+1);
      m_ghostRight(gRow, 10) = m_state(ind+2);
      m_ghostRight(gRow, 11) = m_state(ind+3);
    }

    if (back2 == -1){
      const auto ind = front1*numDofPerCell;
      m_ghostBack(gRow, 8)  = m_state(ind);
      m_ghostBack(gRow, 9)  = m_state(ind+1);
      m_ghostBack(gRow, 10) = -m_state(ind+2);
      m_ghostBack(gRow, 11) = m_state(ind+3);
    }
  }

private:
  const int m_stencilSize;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  ghost_t & m_ghostLeft;
  ghost_t & m_ghostFront;
  ghost_t & m_ghostRight;
  ghost_t & m_ghostBack;
};

}}
#endif
