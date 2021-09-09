
#ifndef PRESSIODEMOAPPS_GHOST_FILLER_EE3D_HPP_
#define PRESSIODEMOAPPS_GHOST_FILLER_EE3D_HPP_

namespace pressiodemoapps{ namespace impl{

template<class state_t, class mesh_t, class ghost_t>
class Ghost3dSedov
{

public:
  Ghost3dSedov() = delete;
  Ghost3dSedov(const int stencilSize,
	       const state_t & stateIn,
	       const mesh_t & meshIn,
	       ghost_t & ghostLeft,
	       ghost_t & ghostRight,
	       ghost_t & ghostBack,
	       ghost_t & ghostFront,
	       ghost_t & ghostBottom,
	       ghost_t & ghostTop)
    : m_stencilSize(stencilSize),
      m_state(stateIn),
      m_meshObj(meshIn),
      m_ghostLeft(ghostLeft),
      m_ghostRight(ghostRight),
      m_ghostBack(ghostBack),
      m_ghostFront(ghostFront),
      m_ghostBottom(ghostBottom),
      m_ghostTop(ghostTop)
  {}

  template<int stencilSize, class index_t>
  typename std::enable_if<stencilSize == 3>::type
  operator()(index_t smPt, int gRow)
  {

    constexpr int numDofPerCell = 5;
    const auto & graph = m_meshObj.graph();
    const auto cellGID = graph(smPt, 0);
    const auto uIndex  = cellGID*numDofPerCell;

    const auto left0  = graph(smPt, 1);
    const auto front0 = graph(smPt, 2);
    const auto right0 = graph(smPt, 3);
    const auto back0  = graph(smPt, 4);
    const auto bot0   = graph(smPt, 5);
    const auto top0   = graph(smPt, 6);

    if (left0 == -1)
    {
      m_ghostLeft(gRow, 0) = m_state(uIndex);
      m_ghostLeft(gRow, 1) = -m_state(uIndex+1);
      m_ghostLeft(gRow, 2) = m_state(uIndex+2);
      m_ghostLeft(gRow, 3) = m_state(uIndex+3);
      m_ghostLeft(gRow, 4) = m_state(uIndex+4);
    }

    if (front0 == -1)
    {
      m_ghostFront(gRow, 0) = m_state(uIndex);
      m_ghostFront(gRow, 1) = m_state(uIndex+1);
      m_ghostFront(gRow, 2) = m_state(uIndex+2);
      m_ghostFront(gRow, 3) = m_state(uIndex+3);
      m_ghostFront(gRow, 4) = m_state(uIndex+4);
    }

    if (right0 == -1)
    {
      m_ghostRight(gRow, 0) = m_state(uIndex);
      m_ghostRight(gRow, 1) = m_state(uIndex+1);
      m_ghostRight(gRow, 2) = m_state(uIndex+2);
      m_ghostRight(gRow, 3) = m_state(uIndex+3);
      m_ghostRight(gRow, 4) = m_state(uIndex+4);
    }

    if (back0 == -1)
    {
      m_ghostBack(gRow, 0) = m_state(uIndex);
      m_ghostBack(gRow, 1) = m_state(uIndex+1);
      m_ghostBack(gRow, 2) = -m_state(uIndex+2);
      m_ghostBack(gRow, 3) = m_state(uIndex+3);
      m_ghostBack(gRow, 4) = m_state(uIndex+4);
    }

    if (top0 == -1)
    {
      m_ghostTop(gRow, 0) = m_state(uIndex);
      m_ghostTop(gRow, 1) = m_state(uIndex+1);
      m_ghostTop(gRow, 2) = m_state(uIndex+2);
      m_ghostTop(gRow, 3) = m_state(uIndex+3);
      m_ghostTop(gRow, 4) = m_state(uIndex+4);
    }

    if (bot0 == -1)
    {
      m_ghostBottom(gRow, 0) = m_state(uIndex);
      m_ghostBottom(gRow, 1) = m_state(uIndex+1);
      m_ghostBottom(gRow, 2) = m_state(uIndex+2);
      m_ghostBottom(gRow, 3) = -m_state(uIndex+3);
      m_ghostBottom(gRow, 4) = m_state(uIndex+4);
    }
  }

  template<int stencilSize, class index_t>
  typename std::enable_if<stencilSize == 5>::type
  operator()(index_t smPt, int gRow)
  {

    constexpr int numDofPerCell = 5;
    const auto & graph = m_meshObj.graph();

    const auto left0   = graph(smPt, 1);
    const auto front0  = graph(smPt, 2);
    const auto right0  = graph(smPt, 3);
    const auto back0   = graph(smPt, 4);
    const auto bot0    = graph(smPt, 5);
    const auto top0    = graph(smPt, 6);
    const auto left1   = graph(smPt, 7);
    const auto front1  = graph(smPt, 8);
    const auto right1  = graph(smPt, 9);
    const auto back1   = graph(smPt, 10);
    const auto bot1    = graph(smPt, 11);
    const auto top1    = graph(smPt, 12);

    const auto left0i  = left0*numDofPerCell;
    const auto front0i = front0*numDofPerCell;
    const auto right0i = right0*numDofPerCell;
    const auto back0i  = back0*numDofPerCell;
    const auto bot0i   = bot0*numDofPerCell;
    const auto top0i   = top0*numDofPerCell;

    this->template operator()<3, index_t>(smPt, gRow);

    // degree 1
    if (left1 == -1)
    {
      m_ghostLeft(gRow, 5) =  m_state(right0i);
      m_ghostLeft(gRow, 6) =  -m_state(right0i+1);
      m_ghostLeft(gRow, 7) =  m_state(right0i+2);
      m_ghostLeft(gRow, 8) =  m_state(right0i+3);
      m_ghostLeft(gRow, 9) =  m_state(right0i+4);
    }

    if (front1 == -1)
    {
      m_ghostFront(gRow, 5) = m_state(back0i);
      m_ghostFront(gRow, 6) = m_state(back0i+1);
      m_ghostFront(gRow, 7) = m_state(back0i+2);
      m_ghostFront(gRow, 8) = m_state(back0i+3);
      m_ghostFront(gRow, 9) = m_state(back0i+4);
    }

    if (right1 == -1)
    {
      m_ghostRight(gRow, 5) = m_state(left0i);
      m_ghostRight(gRow, 6) = m_state(left0i+1);
      m_ghostRight(gRow, 7) = m_state(left0i+2);
      m_ghostRight(gRow, 8) = m_state(left0i+3);
      m_ghostRight(gRow, 9) = m_state(left0i+4);
    }

    if (back1 == -1)
    {
      m_ghostBack(gRow, 5) =  m_state(front0i);
      m_ghostBack(gRow, 6) =  m_state(front0i+1);
      m_ghostBack(gRow, 7) =  -m_state(front0i+2);
      m_ghostBack(gRow, 8) =  m_state(front0i+3);
      m_ghostBack(gRow, 9) =  m_state(front0i+4);
    }

    if (top1 == -1)
    {
      m_ghostTop(gRow, 5) = m_state(bot0i);
      m_ghostTop(gRow, 6) = m_state(bot0i+1);
      m_ghostTop(gRow, 7) = m_state(bot0i+2);
      m_ghostTop(gRow, 8) = m_state(bot0i+3);
      m_ghostTop(gRow, 9) = m_state(bot0i+4);
    }

    if (bot1 == -1)
    {
      m_ghostBottom(gRow, 5) =  m_state(top0i);
      m_ghostBottom(gRow, 6) =  m_state(top0i+1);
      m_ghostBottom(gRow, 7) =  m_state(top0i+2);
      m_ghostBottom(gRow, 8) =  -m_state(top0i+3);
      m_ghostBottom(gRow, 9) =  m_state(top0i+4);
    }
  }

  template<int stencilSize, class index_t>
  typename std::enable_if<stencilSize==7>::type
  operator()(index_t smPt, int gRow)
  {
    throw std::runtime_error("For 3d, only 1st order is currently working");
  }

private:
  const int m_stencilSize;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  ghost_t & m_ghostLeft;
  ghost_t & m_ghostRight;
  ghost_t & m_ghostBack;
  ghost_t & m_ghostFront;
  ghost_t & m_ghostBottom;
  ghost_t & m_ghostTop;
};

}}
#endif




  // template<int stencilSize, class index_t>
  // typename std::enable_if<stencilSize == 7>::type
  // operator()(index_t smPt, int gRow)
  // {
  //   constexpr int numDofPerCell = 5;
  //   const auto & graph = m_meshObj.graph();
  //   const auto cellGID = graph(smPt, 0);
  //   const auto uIndex  = cellGID*numDofPerCell;

  //   this->template operator()<3, index_t>(smPt, gRow);

  //   const auto left0  = graph(smPt, 1);
  //   const auto front0 = graph(smPt, 2);
  //   const auto right0 = graph(smPt, 3);
  //   const auto back0  = graph(smPt, 4);
  //   const auto bot0   = graph(smPt, 5);
  //   const auto top0   = graph(smPt, 6);

  //   const auto left1  = graph(smPt, 7);
  //   const auto front1 = graph(smPt, 8);
  //   const auto right1 = graph(smPt, 9);
  //   const auto back1  = graph(smPt, 10);
  //   const auto bot1   = graph(smPt, 11);
  //   const auto top1   = graph(smPt, 12);

  //   const auto left2  = graph(smPt, 13);
  //   const auto front2 = graph(smPt, 14);
  //   const auto right2 = graph(smPt, 15);
  //   const auto back2  = graph(smPt, 16);
  //   const auto bot2   = graph(smPt, 17);
  //   const auto top2   = graph(smPt, 18);

  //   if (left1 == -1){
  //     const auto ind = right0*numDofPerCell;
  //     m_ghostLeft(gRow, 5) =  m_state(ind);
  //     m_ghostLeft(gRow, 6) =  -m_state(ind+1);
  //     m_ghostLeft(gRow, 7) =  m_state(ind+2);
  //     m_ghostLeft(gRow, 8) =  m_state(ind+3);
  //     m_ghostLeft(gRow, 9) =  m_state(ind+4);
  //   }

  //   if (front1 == -1){
  //     const auto ind = back0*numDofPerCell;
  //     m_ghostFront(gRow, 5) = m_state(ind);
  //     m_ghostFront(gRow, 6) = m_state(ind+1);
  //     m_ghostFront(gRow, 7) = m_state(ind+2);
  //     m_ghostFront(gRow, 8) = m_state(ind+3);
  //     m_ghostFront(gRow, 9) = m_state(ind+4);
  //   }

  //   if (right1 == -1){
  //     const auto ind = left0*numDofPerCell;
  //     m_ghostRight(gRow, 5) = m_state(ind);
  //     m_ghostRight(gRow, 6) = m_state(ind+1);
  //     m_ghostRight(gRow, 7) = m_state(ind+2);
  //     m_ghostRight(gRow, 8) = m_state(ind+3);
  //     m_ghostRight(gRow, 9) = m_state(ind+4);
  //   }

  //   if (back1 == -1){
  //     const auto ind = front0*numDofPerCell;
  //     m_ghostBack(gRow, 5) =  m_state(ind);
  //     m_ghostBack(gRow, 6) =  m_state(ind+1);
  //     m_ghostBack(gRow, 7) =  m_state(ind+2);
  //     m_ghostBack(gRow, 8) =  m_state(ind+3);
  //     m_ghostBack(gRow, 9) =  m_state(ind+4);
  //   }

  //   if (top1 == -1){
  //     const auto ind = bot0*numDofPerCell;
  //     m_ghostTop(gRow, 5) = m_state(ind);
  //     m_ghostTop(gRow, 6) = m_state(ind+1);
  //     m_ghostTop(gRow, 7) = m_state(ind+2);
  //     m_ghostTop(gRow, 8) = m_state(ind+3);
  //     m_ghostTop(gRow, 9) = m_state(ind+4);
  //   }

  //   if (bot1 == -1){
  //     const auto ind = top0*numDofPerCell;
  //     m_ghostBottom(gRow, 5) =  m_state(ind);
  //     m_ghostBottom(gRow, 6) =  m_state(ind+1);
  //     m_ghostBottom(gRow, 7) =  m_state(ind+2);
  //     m_ghostBottom(gRow, 8) =  m_state(ind+3);
  //     m_ghostBottom(gRow, 9) =  m_state(ind+4);
  //   }

  //   // ----
  //   if (left2 == -1){
  //     const auto ind = right1*numDofPerCell;
  //     m_ghostLeft(gRow, 10) =  m_state(ind);
  //     m_ghostLeft(gRow, 11) =  -m_state(ind+1);
  //     m_ghostLeft(gRow, 12) =  m_state(ind+2);
  //     m_ghostLeft(gRow, 13) =  m_state(ind+3);
  //     m_ghostLeft(gRow, 14) =  m_state(ind+4);
  //   }

  //   if (front2 == -1){
  //     const auto ind = back1*numDofPerCell;
  //     m_ghostFront(gRow, 10) = m_state(ind);
  //     m_ghostFront(gRow, 11) = m_state(ind+1);
  //     m_ghostFront(gRow, 12) = m_state(ind+2);
  //     m_ghostFront(gRow, 13) = m_state(ind+3);
  //     m_ghostFront(gRow, 14) = m_state(ind+4);
  //   }

  //   if (right2 == -1){
  //     const auto ind = left1*numDofPerCell;
  //     m_ghostRight(gRow, 10) = m_state(ind);
  //     m_ghostRight(gRow, 11) = m_state(ind+1);
  //     m_ghostRight(gRow, 12) = m_state(ind+2);
  //     m_ghostRight(gRow, 13) = m_state(ind+3);
  //     m_ghostRight(gRow, 14) = m_state(ind+4);
  //   }

  //   if (back2 == -1){
  //     const auto ind = front1*numDofPerCell;
  //     m_ghostBack(gRow, 10) =  m_state(ind);
  //     m_ghostBack(gRow, 11) =  m_state(ind+1);
  //     m_ghostBack(gRow, 12) =  m_state(ind+2);
  //     m_ghostBack(gRow, 13) =  m_state(ind+3);
  //     m_ghostBack(gRow, 14) =  m_state(ind+4);
  //   }

  //   if (top2 == -1){
  //     const auto ind = bot1*numDofPerCell;
  //     m_ghostTop(gRow, 10) = m_state(ind);
  //     m_ghostTop(gRow, 11) = m_state(ind+1);
  //     m_ghostTop(gRow, 12) = m_state(ind+2);
  //     m_ghostTop(gRow, 13) = m_state(ind+3);
  //     m_ghostTop(gRow, 14) = m_state(ind+4);
  //   }

  //   if (bot2 == -1){
  //     const auto ind = top1*numDofPerCell;
  //     m_ghostBottom(gRow, 10) =  m_state(ind);
  //     m_ghostBottom(gRow, 11) =  m_state(ind+1);
  //     m_ghostBottom(gRow, 12) =  m_state(ind+2);
  //     m_ghostBottom(gRow, 13) =  m_state(ind+3);
  //     m_ghostBottom(gRow, 14) =  m_state(ind+4);
  //   }
//}
