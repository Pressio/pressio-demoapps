
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

  template<int stencilSize, class index_t>
  typename std::enable_if<stencilSize == 3>::type
  operator()(index_t smPt, int gRow)
  {

    constexpr int numDofPerCell = 4;
    const auto & graph = m_meshObj.graph();
    const auto cellGID = graph(smPt, 0);
    const auto uIndex  = cellGID*numDofPerCell;

    const auto w0 = graph(smPt, 1);
    const auto n0 = graph(smPt, 2);
    const auto e0 = graph(smPt, 3);
    const auto s0 = graph(smPt, 4);

    if (w0 == -1)
    {
      m_ghostLeft(gRow, 0) = m_state(uIndex);
      m_ghostLeft(gRow, 1) = m_state(uIndex+1);
      m_ghostLeft(gRow, 2) = m_state(uIndex+2);
      m_ghostLeft(gRow, 3) = m_state(uIndex+3);
    }

    if (n0 == -1)
    {
      m_ghostTop(gRow, 0) = m_state(uIndex);
      m_ghostTop(gRow, 1) = m_state(uIndex+1);
      m_ghostTop(gRow, 2) = m_state(uIndex+2);
      m_ghostTop(gRow, 3) = m_state(uIndex+3);
    }

    if (e0 == -1)
    {
      m_ghostRight(gRow, 0) = m_state(uIndex);
      m_ghostRight(gRow, 1) = m_state(uIndex+1);
      m_ghostRight(gRow, 2) = m_state(uIndex+2);
      m_ghostRight(gRow, 3) = m_state(uIndex+3);
    }

    if (s0 == -1)
    {
      m_ghostBottom(gRow, 0) = m_state(uIndex);
      m_ghostBottom(gRow, 1) = m_state(uIndex+1);
      m_ghostBottom(gRow, 2) = m_state(uIndex+2);
      m_ghostBottom(gRow, 3) = m_state(uIndex+3);
    }
  }

  template<int stencilSize, class index_t>
  typename std::enable_if<stencilSize == 7>::type
  operator()(index_t smPt, int gRow)
  {
    constexpr int numDofPerCell = 4;
    const auto & graph = m_meshObj.graph();
    const auto cellGID = graph(smPt, 0);
    const auto uIndex  = cellGID*numDofPerCell;

    this->template operator()<3, index_t>(smPt, gRow);
    const auto w0 = graph(smPt, 1);
    const auto n0 = graph(smPt, 2);
    const auto e0 = graph(smPt, 3);
    const auto s0 = graph(smPt, 4);
    const auto w1 = graph(smPt, 5);
    const auto n1 = graph(smPt, 6);
    const auto e1 = graph(smPt, 7);
    const auto s1 = graph(smPt, 8);
    const auto w2 = graph(smPt, 9);
    const auto n2 = graph(smPt, 10);
    const auto e2 = graph(smPt, 11);
    const auto s2 = graph(smPt, 12);

    if (w1 == -1){
      const auto ind = e0*numDofPerCell;
      m_ghostLeft(gRow, 4) = m_state(ind);
      m_ghostLeft(gRow, 5) = m_state(ind+1);
      m_ghostLeft(gRow, 6) = m_state(ind+2);
      m_ghostLeft(gRow, 7) = m_state(ind+3);
    }

    if (n1 == -1){
      const auto ind = s0*numDofPerCell;
      m_ghostTop(gRow, 4) = m_state(ind);
      m_ghostTop(gRow, 5) = m_state(ind+1);
      m_ghostTop(gRow, 6) = m_state(ind+2);
      m_ghostTop(gRow, 7) = m_state(ind+3);
    }

    if (e1 == -1){
      const auto ind = w0*numDofPerCell;
      m_ghostRight(gRow, 4) = m_state(ind);
      m_ghostRight(gRow, 5) = m_state(ind+1);
      m_ghostRight(gRow, 6) = m_state(ind+2);
      m_ghostRight(gRow, 7) = m_state(ind+3);
    }

    if (s1 == -1){
      const auto ind = n0*numDofPerCell;
      m_ghostBottom(gRow, 4) = m_state(ind);
      m_ghostBottom(gRow, 5) = m_state(ind+1);
      m_ghostBottom(gRow, 6) = m_state(ind+2);
      m_ghostBottom(gRow, 7) = m_state(ind+3);
    }

    if (w2 == -1){
      const auto ind = e1*numDofPerCell;
      m_ghostLeft(gRow, 8)  = m_state(ind);
      m_ghostLeft(gRow, 9)  = m_state(ind+1);
      m_ghostLeft(gRow, 10) = m_state(ind+2);
      m_ghostLeft(gRow, 11) = m_state(ind+3);
    }

    if (n2 == -1){
      const auto ind = s1*numDofPerCell;
      m_ghostTop(gRow, 8) = m_state(ind);
      m_ghostTop(gRow, 9) = m_state(ind+1);
      m_ghostTop(gRow, 10) = m_state(ind+2);
      m_ghostTop(gRow, 11) = m_state(ind+3);
    }

    if (e2 == -1){
      const auto ind = w1*numDofPerCell;
      m_ghostRight(gRow, 8) = m_state(ind);
      m_ghostRight(gRow, 9) = m_state(ind+1);
      m_ghostRight(gRow, 10) = m_state(ind+2);
      m_ghostRight(gRow, 11) = m_state(ind+3);
    }

    if (s2 == -1){
      const auto ind = n1*numDofPerCell;
      m_ghostBottom(gRow, 8) = m_state(ind);
      m_ghostBottom(gRow, 9) = m_state(ind+1);
      m_ghostBottom(gRow, 10) = m_state(ind+2);
      m_ghostBottom(gRow, 11) = m_state(ind+3);
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
