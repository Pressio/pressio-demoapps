
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

  void operator()()
  {
    const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    const auto & graph = m_meshObj.graph();

    // we assume:
    // smPt=0 is at left bd
    // smPt=sampleMeshSize-1 is at right bd
    const auto cellGIDL    = graph(0, 0);
    const auto cellGIDLe0  = graph(0, 2);
    const auto uLe0i	   = cellGIDLe0*3;
    // const auto myXL = x(cellGIDL);
    if (graph(0,1) != -1){
      throw std::runtime_error("Point not on boundary, something wrong");
    }

    const auto cellGIDR   = graph(sampleMeshSize-1, 0);
    const auto cellGIDRw0 = graph(sampleMeshSize-1, 1);
    const auto uRw0i	  = cellGIDRw0*3;
    // const auto myXR = x(cellGIDR);
    if (graph(sampleMeshSize-1, 2) != -1){
      throw std::runtime_error("Point not on boundary, something wrong");
    }

    // this is the 0-degree ghost cell
    m_ghostLeft(0)  = m_state(cellGIDL*3);
    m_ghostLeft(1)  = m_state(cellGIDL*3+1);
    m_ghostLeft(2)  = m_state(cellGIDL*3+2);

    m_ghostRight(0) = m_state(cellGIDR*3);
    m_ghostRight(1) = m_state(cellGIDR*3+1);
    m_ghostRight(2) = m_state(cellGIDR*3+2);

    if (m_stencilSize==5)
      {
	m_ghostLeft(3) = m_state(uLe0i);
	m_ghostLeft(4) = m_state(uLe0i+1);
	m_ghostLeft(5) = m_state(uLe0i+2);

	m_ghostRight(3) = m_state(uRw0i);
	m_ghostRight(4) = m_state(uRw0i+1);
	m_ghostRight(5) = m_state(uRw0i+2);
      }

    if (m_stencilSize==7)
      {
	const auto cellGIDLe1 = graph(0, 4);
	const auto cellGIDRw1 = graph(sampleMeshSize-1, 3);
	m_ghostLeft(6) = m_state(cellGIDLe1*3);
	m_ghostLeft(7) = m_state(cellGIDLe1*3+1);
	m_ghostLeft(8) = m_state(cellGIDLe1*3+2);
	m_ghostLeft(3) = m_state(uLe0i);
	m_ghostLeft(4) = m_state(uLe0i+1);
	m_ghostLeft(5) = m_state(uLe0i+2);

	m_ghostRight(3) = m_state(uRw0i);
	m_ghostRight(4) = m_state(uRw0i+1);
	m_ghostRight(5) = m_state(uRw0i+2);
	m_ghostRight(6) = m_state(cellGIDRw1*3);
	m_ghostRight(7) = m_state(cellGIDRw1*3+1);
	m_ghostRight(8) = m_state(cellGIDRw1*3+2);
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
