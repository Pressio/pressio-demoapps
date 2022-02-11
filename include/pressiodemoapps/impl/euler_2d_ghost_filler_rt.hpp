
#ifndef PRESSIODEMOAPPS_GHOST_FILLER_EE2D_RAYTAY_INSTABILITY_HPP_
#define PRESSIODEMOAPPS_GHOST_FILLER_EE2D_RAYTAY_INSTABILITY_HPP_

namespace pressiodemoapps{ namespace impl{

template<class state_t, class mesh_t, class ghost_t>
class RayTay2dGhostFiller
{

  using scalar_type = typename mesh_t::scalar_t;

public:
  RayTay2dGhostFiller() = delete;
  RayTay2dGhostFiller(const int stencilSize,
		      const state_t & stateIn,
		      const mesh_t & meshIn,
		      const scalar_type gamma,
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
  {
    constexpr auto zero = static_cast<scalar_type>(0);
    constexpr auto one  = static_cast<scalar_type>(1);
    constexpr auto two  = static_cast<scalar_type>(2);
    constexpr auto five = static_cast<scalar_type>(5);
    const auto gammaMinusOneInv = one/(gamma-one);

    std::array<scalar_type, 4> prim;

    // compute conservative BCs at bottom
    prim[0] = two;
    prim[1] = zero;
    prim[2] = zero;
    prim[3] = one;
    conservBot[0] = prim[0];
    conservBot[1] = prim[0]*prim[1];
    conservBot[2] = prim[0]*prim[2];
    conservBot[3] = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv,
							      prim);
    // compute conservative BCs at top
    prim[0] = one;
    prim[1] = zero;
    prim[2] = zero;
    prim[3] = five/two;
    conservTop[0] = prim[0];
    conservTop[1] = prim[0]*prim[1];
    conservTop[2] = prim[0]*prim[2];
    conservTop[3] = eulerEquationsComputeEnergyFromPrimitive2(gammaMinusOneInv,
							      prim);
  }

  template<class index_t>
  void operator()(index_t smPt, int gRow)
  {
    constexpr int numDofPerCell = 4;

    const auto & graph = m_meshObj.graph();
    const auto cellGID = graph(smPt, 0);
    const auto uIndex  = cellGID*numDofPerCell;
    const auto left0   = graph(smPt, 1);
    const auto front0  = graph(smPt, 2);
    const auto right0  = graph(smPt, 3);
    const auto back0   = graph(smPt, 4);

    if (left0 == -1)
    {
      m_ghostLeft(gRow, 0) = m_state(uIndex+0);
      m_ghostLeft(gRow, 1) = -m_state(uIndex+1);
      m_ghostLeft(gRow, 2) = m_state(uIndex+2);
      m_ghostLeft(gRow, 3) = m_state(uIndex+3);
    }

    if (front0 == -1){
      m_ghostFront(gRow, 0) = conservTop[0];
      m_ghostFront(gRow, 1) = conservTop[1];
      m_ghostFront(gRow, 2) = conservTop[2];
      m_ghostFront(gRow, 3) = conservTop[3];
    }

    if (right0 == -1)
    {
      m_ghostRight(gRow, 0) = m_state(uIndex);
      m_ghostRight(gRow, 1) = -m_state(uIndex+1);
      m_ghostRight(gRow, 2) = m_state(uIndex+2);
      m_ghostRight(gRow, 3) = m_state(uIndex+3);
    }

    if (back0 == -1)
    {
      m_ghostBack(gRow, 0) = conservBot[0];
      m_ghostBack(gRow, 1) = conservBot[1];
      m_ghostBack(gRow, 2) = conservBot[2];
      m_ghostBack(gRow, 3) = conservBot[3];
    }

    if (m_stencilSize >= 5){
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
	m_ghostFront(gRow, 4) = conservTop[0];
	m_ghostFront(gRow, 5) = conservTop[1];
	m_ghostFront(gRow, 6) = conservTop[2];
	m_ghostFront(gRow, 7) = conservTop[3];
      }

      if (right1 == -1){
	const auto ind = left0*numDofPerCell;
	m_ghostRight(gRow, 4) = m_state(ind);
	m_ghostRight(gRow, 5) = -m_state(ind+1);
	m_ghostRight(gRow, 6) = m_state(ind+2);
	m_ghostRight(gRow, 7) = m_state(ind+3);
      }

      if (back1 == -1){
	const auto ind = front0*numDofPerCell;
	m_ghostBack(gRow, 4) = conservBot[0];
	m_ghostBack(gRow, 5) = conservBot[1];
	m_ghostBack(gRow, 6) = conservBot[2];
	m_ghostBack(gRow, 7) = conservBot[3];
      }
    }

    if (m_stencilSize == 7){
      const auto left1  = graph(smPt, 5);
      const auto front1 = graph(smPt, 6);
      const auto right1 = graph(smPt, 7);
      const auto back1  = graph(smPt, 8);
      const auto left2  = graph(smPt, 9);
      const auto front2 = graph(smPt, 10);
      const auto right2 = graph(smPt, 11);
      const auto back2  = graph(smPt, 12);

      if (left2 == -1){
	const auto ind = right1*numDofPerCell;
	m_ghostLeft(gRow, 8)  = m_state(ind);
	m_ghostLeft(gRow, 9)  = -m_state(ind+1);
	m_ghostLeft(gRow, 10) = m_state(ind+2);
	m_ghostLeft(gRow, 11) = m_state(ind+3);
      }

      if (front2 == -1){
	const auto ind = back1*numDofPerCell;
	m_ghostFront(gRow, 8)  = conservTop[0];
	m_ghostFront(gRow, 9)  = conservTop[1];
	m_ghostFront(gRow, 10) = conservTop[2];
	m_ghostFront(gRow, 11) = conservTop[3];
      }

      if (right2 == -1){
	const auto ind = left1*numDofPerCell;
	m_ghostRight(gRow, 8)  = m_state(ind);
	m_ghostRight(gRow, 9)  = -m_state(ind+1);
	m_ghostRight(gRow, 10) = m_state(ind+2);
	m_ghostRight(gRow, 11) = m_state(ind+3);
      }

      if (back2 == -1){
	const auto ind = front1*numDofPerCell;
	m_ghostBack(gRow, 8)  = conservBot[0];
	m_ghostBack(gRow, 9)  = conservBot[1];
	m_ghostBack(gRow, 10) = conservBot[2];
	m_ghostBack(gRow, 11) = conservBot[3];
      }
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

  // store conservative Dirichlet BCs top/bottom
  std::array<scalar_type, 4> conservBot;
  std::array<scalar_type, 4> conservTop;
};

}}
#endif
